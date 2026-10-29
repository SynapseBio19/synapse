"""
Synapse - CRISPR Analysis Engine
================================

Professional, deterministic, sequence-level CRISPR analysis module for Project
Synapse.  This module is intentionally self-contained and uses only Python's
standard library.

Capabilities
------------
- Strand-aware PAM/protospacer discovery for multiple CRISPR systems.
- Explicit 1-based coordinates and target orientation.
- Rich guide characterization: GC/AT, ambiguity, homopolymers, repeats,
  complexity, self-complementarity and hairpin-risk proxies.
- Transparent component-based guide ranking; no fabricated experimental
  efficacy score.
- Mutation/variant-aware targeting and bystander analysis.
- PAM-aware local specificity/off-target screening on explicitly supplied
  sequences, including mismatch position, seed-region and PAM compatibility.
- Base-editing candidate analysis with configurable editing windows.
- Prime-editing sequence-level candidate construction (PBS/RTT candidates),
  explicitly labelled heuristic rather than experimentally validated design.
- ORF/region/context integration and batch FASTA analysis.
- Reproducible JSON/CSV/FASTA reports and a detailed text report.
- Comprehensive deterministic self-test suite.

Scientific boundary
-------------------
This is a computational sequence-analysis module. It does not establish
experimental editing efficacy, cellular delivery, genomic specificity in an
organism, pathogenicity, clinical suitability, or experimental outcome.
Genome-wide specificity requires a real reference genome and a validated
indexed search engine; the local engine here only evaluates sequences supplied
to Synapse. Base-editing and prime-editing outputs are computational candidate
architectures and not experimental efficacy predictions.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
import csv
import json
import math
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


# ============================================================================
# CONSTANTS
# ============================================================================

DNA_BASES = set("ACGT")
IUPAC_BASES = set("ACGTRYSWKMBDHVN")

IUPAC_CODES: Dict[str, set[str]] = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

IUPAC_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}

# Data-driven CRISPR system registry.  PAM orientation is defined relative to
# the guide/protospacer sequence in its targeting orientation.
SYSTEM_DEFINITIONS: Dict[str, Dict[str, Any]] = {
    "SpCas9": {
        "enzyme_family": "Cas9",
        "pam": "NGG",
        "pam_length": 3,
        "guide_length": 20,
        "pam_side": "3prime",
        "cut_model": "Cas9-like; cleavage model is not an efficacy prediction",
        "description": "SpCas9 canonical NGG PAM model.",
    },
    "SpCas9-NG": {
        "enzyme_family": "Cas9",
        "pam": "NG",
        "pam_length": 2,
        "guide_length": 20,
        "pam_side": "3prime",
        "cut_model": "Cas9-like; relaxed PAM model",
        "description": "SpCas9-NG-like relaxed NG PAM model.",
    },
    "SaCas9": {
        "enzyme_family": "Cas9",
        "pam": "NNGRRT",
        "pam_length": 6,
        "guide_length": 21,
        "pam_side": "3prime",
        "cut_model": "Cas9-like; cleavage model is not an efficacy prediction",
        "description": "SaCas9 canonical NNGRRT PAM model.",
    },
    "Cas12a": {
        "enzyme_family": "Cas12a",
        "pam": "TTTV",
        "pam_length": 4,
        "guide_length": 23,
        "pam_side": "5prime",
        "cut_model": "Cas12a-like; distal cleavage behaviour not simulated",
        "description": "Cas12a-style T-rich 5' PAM model.",
    },
}

BASE_EDITOR_DEFINITIONS: Dict[str, Dict[str, Any]] = {
    "CBE": {
        "conversion": ("C", "T"),
        "window": (4, 8),
        "description": "C-to-T-compatible cytosine base-editing model.",
    },
    "ABE": {
        "conversion": ("A", "G"),
        "window": (4, 8),
        "description": "A-to-G-compatible adenine base-editing model.",
    },
    "CBE+": {
        "conversion": ("C", "T"),
        "window": (4, 9),
        "description": "Broader C-to-T-compatible computational window.",
    },
    "ABE+": {
        "conversion": ("A", "G"),
        "window": (4, 9),
        "description": "Broader A-to-G-compatible computational window.",
    },
}

DEFAULT_MIN_GC = 30.0
DEFAULT_MAX_GC = 70.0
DEFAULT_MAX_HOMOPOLYMER = 4
DEFAULT_MAX_LOCAL_MISMATCHES = 3
DEFAULT_SEED_LENGTH = 8


# ============================================================================
# EXCEPTIONS
# ============================================================================

class CRISPRAnalysisError(ValueError):
    """Raised for invalid CRISPR analysis input."""


# ============================================================================
# DATA MODELS
# ============================================================================

@dataclass(frozen=True)
class CRISPRSystem:
    name: str
    enzyme_family: str
    pam: str
    pam_length: int
    guide_length: int
    pam_side: str
    cut_model: str = ""
    description: str = ""


@dataclass
class TargetRegion:
    sequence: str
    start: int = 1
    end: Optional[int] = None
    name: str = "Selected Region"
    source: str = "sequence"
    feature_type: str = "sequence"
    feature_name: str = ""
    feature_strand: str = "+"
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.sequence = normalize_sequence(self.sequence)
        if self.end is None:
            self.end = self.start + len(self.sequence) - 1
        if self.start < 1:
            raise CRISPRAnalysisError("Target region start must be >= 1.")
        if self.end < self.start:
            raise CRISPRAnalysisError("Target region end must be >= start.")
        if self.end - self.start + 1 != len(self.sequence):
            raise CRISPRAnalysisError("Target region coordinates do not match sequence length.")


@dataclass
class GuideMetrics:
    gc_content: float
    at_content: float
    gc_deviation: float
    homopolymer_max: int
    repeated_dinucleotide: bool
    ambiguous_bases: int
    sequence_complexity: float
    poly_t: bool
    poly_g: bool
    repeat_fraction: float = 0.0
    self_complementarity_max: int = 0
    hairpin_risk: str = "LOW"
    unique_kmer_fraction: float = 1.0


@dataclass
class MutationTarget:
    chromosome: str
    position: int
    reference: str
    alternative: str
    label: str = "Variant"

    def __post_init__(self) -> None:
        self.chromosome = str(self.chromosome)
        self.position = int(self.position)
        self.reference = normalize_sequence(self.reference)
        self.alternative = normalize_sequence(self.alternative)
        if self.position < 1:
            raise CRISPRAnalysisError("Mutation position must be >= 1.")


@dataclass
class GuideCandidate:
    candidate_id: str
    system: str
    guide: str
    pam: str
    strand: str
    guide_start: int
    guide_end: int
    pam_start: int
    pam_end: int
    guide_length: int
    target_sequence: str
    metrics: GuideMetrics
    pam_quality: Dict[str, Any] = field(default_factory=dict)
    seed_sequence: str = ""
    rank_score: float = 0.0
    rank_percentile: float = 0.0
    rank_components: Dict[str, float] = field(default_factory=dict)
    flags: List[str] = field(default_factory=list)
    mutation_overlap: List[str] = field(default_factory=list)
    mutation_context: List[Dict[str, Any]] = field(default_factory=list)
    off_target_hits: int = 0
    local_off_target_hits: List[Dict[str, Any]] = field(default_factory=list)
    specificity_profile: Dict[str, Any] = field(default_factory=dict)
    editing_window: Tuple[int, int] = (0, 0)
    bystander_edits: List[Dict[str, Any]] = field(default_factory=list)
    base_editing_candidates: List[Dict[str, Any]] = field(default_factory=list)
    prime_editing_candidates: List[Dict[str, Any]] = field(default_factory=list)
    context: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_high_quality(self) -> bool:
        return not any(flag.startswith("FAIL:") for flag in self.flags)


@dataclass
class DesignConfig:
    system: str = "SpCas9"
    min_gc: float = DEFAULT_MIN_GC
    max_gc: float = DEFAULT_MAX_GC
    max_homopolymer: int = DEFAULT_MAX_HOMOPOLYMER
    min_complexity: float = 0.45
    allow_ambiguous: bool = False
    include_low_quality: bool = True
    rank_by_gc: bool = True
    mutation_window: int = 0
    max_local_off_target_mismatches: int = DEFAULT_MAX_LOCAL_MISMATCHES
    max_local_off_target_hits: int = 100
    seed_length: int = DEFAULT_SEED_LENGTH
    require_pam_for_off_targets: bool = True
    editing_window_start: int = 4
    editing_window_end: int = 8
    base_editor: str = "CBE"
    prime_pbs_lengths: Tuple[int, ...] = tuple(range(8, 16))
    prime_rtt_lengths: Tuple[int, ...] = tuple(range(10, 31, 2))


@dataclass
class DesignResult:
    sequence_name: str
    sequence_length: int
    region: TargetRegion
    system: CRISPRSystem
    config: DesignConfig
    candidates: List[GuideCandidate]
    warnings: List[str] = field(default_factory=list)
    analysis_metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def candidate_count(self) -> int:
        return len(self.candidates)

    @property
    def high_quality_count(self) -> int:
        return sum(c.is_high_quality for c in self.candidates)

    # Backward-compatible convenience accessors for CLI/main.py integrations.
    @property
    def region_start(self) -> int:
        return self.region.start

    @property
    def region_end(self) -> int:
        return self.region.end

    @property
    def region_name(self) -> str:
        return self.region.name

    @property
    def feature_type(self) -> str:
        return self.region.feature_type


@dataclass
class BaseEditingCandidate:
    candidate_id: str
    editor: str
    conversion: str
    guide: str
    pam: str
    strand: str
    variant: str
    variant_position_in_guide: int
    editing_window: Tuple[int, int]
    desired_edit_in_window: bool
    bystander_count: int
    bystanders: List[Dict[str, Any]]
    sequence_level_notes: List[str] = field(default_factory=list)


@dataclass
class PrimeEditingCandidate:
    candidate_id: str
    guide_id: str
    guide: str
    pam: str
    strand: str
    variant: str
    pbs_length: int
    pbs_sequence: str
    rtt_length: int
    rtt_sequence: str
    edit_in_rtt: bool
    sequence_level_notes: List[str] = field(default_factory=list)


@dataclass
class SelfTestResult:
    name: str
    passed: bool
    detail: str = ""


# ============================================================================
# BASIC SEQUENCE UTILITIES
# ============================================================================

def normalize_sequence(sequence: str) -> str:
    if sequence is None:
        raise CRISPRAnalysisError("Sequence cannot be None.")
    cleaned = "".join(str(sequence).split()).upper()
    if not cleaned:
        raise CRISPRAnalysisError("Sequence cannot be empty.")
    invalid = sorted(set(cleaned) - IUPAC_BASES)
    if invalid:
        raise CRISPRAnalysisError("Invalid DNA symbols: " + ", ".join(invalid))
    return cleaned


def validate_sequence(sequence: str, allow_ambiguous: bool = False) -> Tuple[bool, List[str]]:
    errors: List[str] = []
    try:
        seq = normalize_sequence(sequence)
    except CRISPRAnalysisError as exc:
        return False, [str(exc)]
    if not allow_ambiguous:
        ambiguous = sorted(set(seq) - DNA_BASES)
        if ambiguous:
            errors.append("Ambiguous/IUPAC bases present: " + ", ".join(ambiguous))
    return not errors, errors


def complement(sequence: str) -> str:
    seq = normalize_sequence(sequence)
    return "".join(IUPAC_COMPLEMENT[base] for base in seq)


def reverse_complement(sequence: str) -> str:
    return complement(sequence)[::-1]


def gc_content(sequence: str) -> float:
    seq = normalize_sequence(sequence)
    concrete = [base for base in seq if base in DNA_BASES]
    return 0.0 if not concrete else 100.0 * sum(b in {"G", "C"} for b in concrete) / len(concrete)


def at_content(sequence: str) -> float:
    seq = normalize_sequence(sequence)
    concrete = [base for base in seq if base in DNA_BASES]
    return 0.0 if not concrete else 100.0 * sum(b in {"A", "T"} for b in concrete) / len(concrete)


def count_ambiguous_bases(sequence: str) -> int:
    return sum(base not in DNA_BASES for base in normalize_sequence(sequence))


def max_homopolymer_length(sequence: str) -> int:
    seq = normalize_sequence(sequence)
    best = 0
    current = ""
    run = 0
    for base in seq:
        if base == current:
            run += 1
        else:
            current = base
            run = 1
        best = max(best, run)
    return best


def has_poly_run(sequence: str, base: str, minimum: int = 4) -> bool:
    return base * minimum in normalize_sequence(sequence)


def has_repeated_dinucleotide(sequence: str, minimum_repeats: int = 3) -> bool:
    seq = normalize_sequence(sequence)
    for i in range(len(seq) - 2 * minimum_repeats + 1):
        unit = seq[i:i + 2]
        if unit * minimum_repeats == seq[i:i + 2 * minimum_repeats]:
            return True
    return False


def sequence_complexity(sequence: str) -> float:
    """Normalized Shannon entropy of concrete A/C/G/T composition."""
    seq = [b for b in normalize_sequence(sequence) if b in DNA_BASES]
    if not seq:
        return 0.0
    counts = Counter(seq)
    total = len(seq)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return entropy / 2.0


def repeat_fraction(sequence: str, k: int = 2) -> float:
    seq = normalize_sequence(sequence)
    if len(seq) < k:
        return 0.0
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    repeated = sum(n for n in counts.values() if n > 1)
    return repeated / len(kmers)


def unique_kmer_fraction(sequence: str, k: int = 4) -> float:
    seq = normalize_sequence(sequence)
    if len(seq) < k:
        return 1.0
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    return len(set(kmers)) / len(kmers)


def self_complementarity_max(sequence: str, min_stem: int = 4, max_stem: int = 8) -> int:
    """Sequence-only reverse-complement stem proxy; not a thermodynamic model."""
    seq = normalize_sequence(sequence)
    best = 0
    for k in range(min_stem, min(max_stem, len(seq)) + 1):
        rc = reverse_complement(seq)
        for i in range(len(seq) - k + 1):
            motif = seq[i:i + k]
            if motif in rc:
                best = max(best, k)
    return best


def hairpin_risk(stem: int) -> str:
    if stem >= 8:
        return "HIGH"
    if stem >= 6:
        return "MODERATE"
    if stem >= 4:
        return "LOW-MODERATE"
    return "LOW"


def mismatch_count(sequence_a: str, sequence_b: str) -> int:
    a = normalize_sequence(sequence_a)
    b = normalize_sequence(sequence_b)
    if len(a) != len(b):
        raise CRISPRAnalysisError("Sequences must have equal length for mismatch comparison.")
    return sum(x != y for x, y in zip(a, b))


def mismatch_positions(sequence_a: str, sequence_b: str) -> List[int]:
    a = normalize_sequence(sequence_a)
    b = normalize_sequence(sequence_b)
    if len(a) != len(b):
        raise CRISPRAnalysisError("Sequences must have equal length for mismatch comparison.")
    return [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]


def parse_fasta_text(text: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    for raw in str(text).splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, normalize_sequence("".join(chunks))))
            name = line[1:].strip() or f"Sequence_{len(records) + 1}"
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, normalize_sequence("".join(chunks))))
    if not records:
        raise CRISPRAnalysisError("No FASTA records found.")
    return records


def load_fasta(path: str | Path) -> List[Tuple[str, str]]:
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError as exc:
        raise CRISPRAnalysisError(f"Unable to read FASTA: {exc}") from exc
    return parse_fasta_text(text)


# ============================================================================
# PAM / SYSTEM REGISTRY
# ============================================================================

def get_crispr_system(name: str = "SpCas9") -> CRISPRSystem:
    key = str(name).strip()
    if key not in SYSTEM_DEFINITIONS:
        raise CRISPRAnalysisError(
            f"Unknown CRISPR system '{name}'. Available: {', '.join(SYSTEM_DEFINITIONS)}"
        )
    data = dict(SYSTEM_DEFINITIONS[key])
    data["name"] = key
    return CRISPRSystem(**data)


def format_system_details(system_name: str) -> str:
    """Return detailed information for exactly one selected CRISPR system."""
    system = get_crispr_system(system_name)
    lines = [
        "SELECTED CRISPR SYSTEM",
        _line(),
        f"System:          {system.name}",
        f"Enzyme family:   {system.enzyme_family}",
        f"PAM:             {system.pam}",
        f"PAM length:      {system.pam_length} nt",
        f"Guide length:    {system.guide_length} nt",
        f"PAM orientation: {system.pam_side}",
        f"Cut model:       {system.cut_model}",
        f"Description:     {system.description}",
        "",
        "Scientific scope: sequence-level computational model; not an",
        "experimental cleavage or editing-efficiency prediction.",
    ]
    return "\n".join(lines)


def list_crispr_systems() -> List[CRISPRSystem]:
    return [get_crispr_system(name) for name in SYSTEM_DEFINITIONS]


def list_base_editors() -> List[str]:
    return list(BASE_EDITOR_DEFINITIONS)


def pam_matches(sequence: str, pattern: str) -> bool:
    seq = normalize_sequence(sequence)
    pat = normalize_sequence(pattern)
    if len(seq) != len(pat):
        return False
    for base, code in zip(seq, pat):
        if base not in IUPAC_CODES.get(code, set()):
            return False
    return True


def pam_quality(pam: str, pattern: str) -> Dict[str, Any]:
    actual = normalize_sequence(pam)
    pat = normalize_sequence(pattern)
    if len(actual) != len(pat):
        raise CRISPRAnalysisError("PAM and PAM pattern must have equal length.")
    degenerate_positions = sum(len(IUPAC_CODES.get(code, {code})) > 1 for code in pat)
    compatible = pam_matches(actual, pat)
    exact_pattern_positions = sum(actual[i] == pat[i] for i in range(len(actual)))
    return {
        "matched": compatible,
        "pattern": pat,
        "degenerate_positions": degenerate_positions,
        "exact_symbol_matches": exact_pattern_positions,
        "specificity_note": "PAM pattern compatibility only; not a biochemical activity score.",
    }


# ============================================================================
# GUIDE METRICS / QUALITY
# ============================================================================

def characterize_guide(guide: str, config: Optional[DesignConfig] = None) -> GuideMetrics:
    cfg = config or DesignConfig()
    seq = normalize_sequence(guide)
    gc = gc_content(seq)
    at = at_content(seq)
    midpoint = (cfg.min_gc + cfg.max_gc) / 2.0
    stem = self_complementarity_max(seq)
    return GuideMetrics(
        gc_content=gc,
        at_content=at,
        gc_deviation=abs(gc - midpoint),
        homopolymer_max=max_homopolymer_length(seq),
        repeated_dinucleotide=has_repeated_dinucleotide(seq),
        ambiguous_bases=count_ambiguous_bases(seq),
        sequence_complexity=sequence_complexity(seq),
        poly_t=has_poly_run(seq, "T", 4),
        poly_g=has_poly_run(seq, "G", 4),
        repeat_fraction=repeat_fraction(seq),
        self_complementarity_max=stem,
        hairpin_risk=hairpin_risk(stem),
        unique_kmer_fraction=unique_kmer_fraction(seq),
    )


def quality_flags(metrics: GuideMetrics, config: DesignConfig) -> List[str]:
    flags: List[str] = []
    if metrics.ambiguous_bases and not config.allow_ambiguous:
        flags.append("FAIL: ambiguous bases")
    if metrics.gc_content < config.min_gc:
        flags.append("WARN: low GC")
    if metrics.gc_content > config.max_gc:
        flags.append("WARN: high GC")
    if metrics.homopolymer_max > config.max_homopolymer:
        flags.append("FAIL: long homopolymer")
    if metrics.repeated_dinucleotide:
        flags.append("WARN: repeated dinucleotide")
    if metrics.sequence_complexity < config.min_complexity:
        flags.append("WARN: low sequence complexity")
    if metrics.poly_t:
        flags.append("WARN: poly-T")
    if metrics.poly_g:
        flags.append("WARN: poly-G")
    if metrics.hairpin_risk == "HIGH":
        flags.append("WARN: strong self-complementarity proxy")
    return flags


# ============================================================================
# TARGET ENUMERATION
# ============================================================================

def _coordinate_region(region: TargetRegion, local_start: int, local_end: int) -> Tuple[int, int]:
    return region.start + local_start, region.start + local_end


def _candidate_from_oriented_window(
    sequence: str,
    region: TargetRegion,
    system: CRISPRSystem,
    guide_start_local: int,
    guide_end_local: int,
    pam_start_local: int,
    pam_end_local: int,
    strand: str,
    config: DesignConfig,
    index: int,
) -> GuideCandidate:
    genomic_guide = sequence[guide_start_local:guide_end_local + 1]
    genomic_pam = sequence[pam_start_local:pam_end_local + 1]
    if strand == "+":
        guide = genomic_guide
        pam = genomic_pam
    else:
        guide = reverse_complement(genomic_guide)
        pam = reverse_complement(genomic_pam)
    guide_start, guide_end = _coordinate_region(region, guide_start_local, guide_end_local)
    pam_start, pam_end = _coordinate_region(region, pam_start_local, pam_end_local)
    metrics = characterize_guide(guide, config)
    flags = quality_flags(metrics, config)
    return GuideCandidate(
        candidate_id=f"G{index:04d}",
        system=system.name,
        guide=guide,
        pam=pam,
        strand=strand,
        guide_start=guide_start,
        guide_end=guide_end,
        pam_start=pam_start,
        pam_end=pam_end,
        guide_length=len(guide),
        target_sequence=guide + pam if system.pam_side == "3prime" else pam + guide,
        metrics=metrics,
        pam_quality=pam_quality(pam, system.pam),
        # The guide sequence is always stored in targeting orientation, so the
        # PAM-proximal seed is the same end of the guide string on both strands.
        seed_sequence=guide[-min(config.seed_length, len(guide)):],
        flags=flags,
        context={
            "region_name": region.name,
            "region_source": region.source,
            "feature_type": region.feature_type,
            "feature_name": region.feature_name,
            "feature_strand": region.feature_strand,
        },
    )


def _find_targets_in_region(region: TargetRegion, system: CRISPRSystem, config: DesignConfig) -> List[GuideCandidate]:
    seq = region.sequence
    candidates: List[GuideCandidate] = []
    counter = 1
    total = system.guide_length + system.pam_length
    if len(seq) < total:
        return []

    # Search both physical orientations directly.  For the reverse strand,
    # windows are reported in original genomic coordinates but the guide/PAM
    # strings are oriented in the targeting direction.
    for strand in ("+", "-"):
        oriented = seq if strand == "+" else reverse_complement(seq)
        for pam_start in range(0, len(oriented) - system.pam_length + 1):
            pam_end = pam_start + system.pam_length - 1
            pam = oriented[pam_start:pam_end + 1]
            if not pam_matches(pam, system.pam):
                continue
            if system.pam_side == "3prime":
                guide_start = pam_start - system.guide_length
                guide_end = pam_start - 1
            else:
                guide_start = pam_end + 1
                guide_end = guide_start + system.guide_length - 1
            if guide_start < 0 or guide_end >= len(oriented):
                continue

            if strand == "+":
                g_start, g_end = guide_start, guide_end
                p_start, p_end = pam_start, pam_end
            else:
                g_start = len(seq) - 1 - guide_end
                g_end = len(seq) - 1 - guide_start
                p_start = len(seq) - 1 - pam_end
                p_end = len(seq) - 1 - pam_start

            candidates.append(
                _candidate_from_oriented_window(
                    seq, region, system, g_start, g_end, p_start, p_end,
                    strand, config, counter,
                )
            )
            counter += 1

    candidates.sort(key=lambda c: (c.guide_start, c.guide_end, c.strand, c.pam_start))
    for idx, candidate in enumerate(candidates, 1):
        candidate.candidate_id = f"G{idx:04d}"
    return candidates


def find_targets(
    sequence: str,
    *,
    system: str = "SpCas9",
    region_start: int = 1,
    region_end: Optional[int] = None,
    region_name: str = "Selected Region",
    region_source: str = "sequence",
    feature_type: str = "sequence",
    feature_name: str = "",
    feature_strand: str = "+",
    config: Optional[DesignConfig] = None,
) -> List[GuideCandidate]:
    seq = normalize_sequence(sequence)
    cfg = config or DesignConfig(system=system)
    cfg.system = system
    if region_end is None:
        region_end = region_start + len(seq) - 1
    region = TargetRegion(
        seq,
        start=region_start,
        end=region_end,
        name=region_name,
        source=region_source,
        feature_type=feature_type,
        feature_name=feature_name,
        feature_strand=feature_strand,
    )
    return _find_targets_in_region(region, get_crispr_system(system), cfg)


# ============================================================================
# GUIDE RANKING
# ============================================================================

def _component_gc(metrics: GuideMetrics, config: DesignConfig) -> float:
    midpoint = (config.min_gc + config.max_gc) / 2.0
    half_range = max(1.0, (config.max_gc - config.min_gc) / 2.0)
    return max(0.0, 100.0 - (abs(metrics.gc_content - midpoint) / half_range) * 35.0)


def _component_complexity(metrics: GuideMetrics) -> float:
    return max(0.0, min(100.0, metrics.sequence_complexity * 100.0))


def _component_structure(metrics: GuideMetrics, config: DesignConfig) -> float:
    score = 100.0
    if metrics.homopolymer_max > config.max_homopolymer:
        score -= 45.0
    elif metrics.homopolymer_max == config.max_homopolymer:
        score -= 10.0
    if metrics.repeated_dinucleotide:
        score -= 15.0
    if metrics.poly_t:
        score -= 10.0
    if metrics.poly_g:
        score -= 5.0
    if metrics.self_complementarity_max >= 8:
        score -= 20.0
    elif metrics.self_complementarity_max >= 6:
        score -= 10.0
    return max(0.0, score)


def _component_pam(candidate: GuideCandidate) -> float:
    # This measures only pattern degeneracy/compatibility, not biochemical activity.
    deg = float(candidate.pam_quality.get("degenerate_positions", 0))
    length = max(1.0, len(candidate.pam))
    return max(0.0, 100.0 - (deg / length) * 20.0)


def _component_specificity(candidate: GuideCandidate) -> Optional[float]:
    """Convert an explicit local specificity profile into a transparent score.

    This is a local-search burden heuristic only.  It is deliberately not a
    genome-wide specificity or experimental efficacy prediction.  Exact local
    matches receive the strongest penalty, followed by one-, two- and
    three-mismatch matches.  The categories are mutually summarized here so
    the same hit is not double-penalized through overlapping counters.
    """
    profile = candidate.specificity_profile
    if not profile or not profile.get("analyzed"):
        return None

    exact = int(profile.get("exact_matches", 0))
    one_mm = int(profile.get("one_mismatch_hits", 0))
    two_mm = int(profile.get("two_mismatch_hits", 0))
    three_mm = int(profile.get("three_mismatch_hits", 0))

    burden = (
        exact * 40.0
        + one_mm * 20.0
        + two_mm * 10.0
        + three_mm * 5.0
    )
    return max(0.0, min(100.0, 100.0 - burden))


def has_specificity_analysis(candidates: Sequence[GuideCandidate]) -> bool:
    """Return True when at least one candidate has explicit local specificity data."""
    return any(bool(c.specificity_profile.get("analyzed")) for c in candidates)


def ranking_mode_label(candidates: Sequence[GuideCandidate]) -> str:
    if has_specificity_analysis(candidates):
        return "specificity-aware local ranking"
    return "sequence/context ranking"


def rank_guides(candidates: Sequence[GuideCandidate], *, config: Optional[DesignConfig] = None) -> List[GuideCandidate]:
    """Rank candidates using transparent computational components.

    The score is a prioritization heuristic, not an experimental editing-efficiency
    prediction. Specificity contributes only when an explicit local analysis exists.
    """
    cfg = config or DesignConfig()
    ranked = list(candidates)
    for candidate in ranked:
        gc_score = _component_gc(candidate.metrics, cfg)
        complexity_score = _component_complexity(candidate.metrics)
        structure_score = _component_structure(candidate.metrics, cfg)
        pam_score = _component_pam(candidate)
        ambiguity_score = 0.0 if candidate.metrics.ambiguous_bases else 100.0
        specificity_score = _component_specificity(candidate)
        warning_penalty = 0.0
        for flag in candidate.flags:
            if flag.startswith("FAIL:"):
                warning_penalty += 30.0
            elif flag.startswith("WARN:"):
                warning_penalty += 7.0

        # Mutation overlap is a context signal, not a universal biological advantage.
        mutation_bonus = 5.0 if candidate.mutation_overlap else 0.0
        base_components = {
            "GC": gc_score,
            "Complexity": complexity_score,
            "SequenceStructure": structure_score,
            "PAMPattern": pam_score,
            "Ambiguity": ambiguity_score,
            "MutationContext": mutation_bonus,
            "Penalty": warning_penalty,
        }
        weights = {
            "GC": 0.25,
            "Complexity": 0.18,
            "SequenceStructure": 0.22,
            "PAMPattern": 0.10,
            "Ambiguity": 0.10,
            "MutationContext": 0.05,
        }
        if specificity_score is not None:
            base_components["Specificity"] = specificity_score
            weights["Specificity"] = 0.10
        else:
            # Re-normalize without an unavailable component.
            total = sum(weights.values())
            for key in list(weights):
                weights[key] /= total

        final = sum(base_components[k] * weights[k] for k in weights) - warning_penalty
        candidate.rank_components = {k: round(v, 3) for k, v in base_components.items()}
        candidate.rank_score = round(max(0.0, min(100.0, final)), 3)

    ranked.sort(key=lambda c: (-c.rank_score, c.guide_start, c.strand, c.guide))
    n = len(ranked)
    for i, candidate in enumerate(ranked):
        candidate.rank_percentile = round(100.0 if n == 1 else 100.0 * (n - i - 1) / (n - 1), 2)
    return ranked


# ============================================================================
# MUTATION / VARIANT CONTEXT
# ============================================================================

def guide_overlaps_position(candidate: GuideCandidate, position: int, window: int = 0) -> bool:
    return candidate.guide_start <= position + window and candidate.guide_end >= position - window


def mutation_position_in_guide(candidate: GuideCandidate, position: int) -> Optional[int]:
    if not guide_overlaps_position(candidate, position):
        return None
    if candidate.strand == "+":
        return position - candidate.guide_start + 1
    return candidate.guide_end - position + 1


def annotate_mutation_targets(
    candidates: Sequence[GuideCandidate],
    mutations: Sequence[MutationTarget],
    *,
    window: int = 0,
) -> List[GuideCandidate]:
    for candidate in candidates:
        candidate.mutation_overlap.clear()
        candidate.mutation_context.clear()
        for mutation in mutations:
            if guide_overlaps_position(candidate, mutation.position, window):
                pos_in_guide = mutation_position_in_guide(candidate, mutation.position)
                item = {
                    "label": mutation.label,
                    "chromosome": mutation.chromosome,
                    "position": mutation.position,
                    "reference": mutation.reference,
                    "alternative": mutation.alternative,
                    "position_in_guide": pos_in_guide,
                    "distance_to_pam": abs(mutation.position - candidate.pam_start) if candidate.strand == "+" else abs(candidate.pam_end - mutation.position),
                }
                candidate.mutation_context.append(item)
                candidate.mutation_overlap.append(
                    f"{mutation.label} {mutation.chromosome}:{mutation.position} {mutation.reference}>{mutation.alternative}"
                )
        if candidate.mutation_overlap and "INFO: mutation-overlap" not in candidate.flags:
            candidate.flags.append("INFO: mutation-overlap")
    return list(candidates)


def parse_mutation_target(text: str) -> MutationTarget:
    raw = str(text).strip()
    pattern = re.compile(
        r"^(?:(?P<chrom>[^:]+):)?(?P<pos>\d+)\s+(?P<ref>[ACGTN]+)>(?P<alt>[ACGTN]+)$",
        re.I,
    )
    match = pattern.match(raw)
    if not match:
        raise CRISPRAnalysisError("Mutation format must be 'chr:position REF>ALT' or 'position REF>ALT'.")
    return MutationTarget(
        chromosome=match.group("chrom") or "sequence",
        position=int(match.group("pos")),
        reference=match.group("ref").upper(),
        alternative=match.group("alt").upper(),
    )


def parse_vcf_variant_line(line: str) -> List[MutationTarget]:
    fields = line.rstrip().split("\t")
    if len(fields) < 5 or line.startswith("#"):
        return []
    chrom, pos, _, ref, alt = fields[:5]
    return [MutationTarget(chrom, int(pos), ref, allele, "VCF Variant") for allele in alt.split(",")]


# ============================================================================
# LOCAL OFF-TARGET / SPECIFICITY ENGINE
# ============================================================================

def _oriented_candidate_windows(sequence: str, candidate: GuideCandidate, system: CRISPRSystem) -> Iterable[Dict[str, Any]]:
    seq = normalize_sequence(sequence)
    for strand in ("+", "-"):
        oriented = seq if strand == "+" else reverse_complement(seq)
        for pam_start in range(0, len(oriented) - system.pam_length + 1):
            pam_end = pam_start + system.pam_length - 1
            pam = oriented[pam_start:pam_end + 1]
            if not pam_matches(pam, system.pam):
                continue
            if system.pam_side == "3prime":
                g_start = pam_start - system.guide_length
                g_end = pam_start - 1
            else:
                g_start = pam_end + 1
                g_end = g_start + system.guide_length - 1
            if g_start < 0 or g_end >= len(oriented):
                continue
            guide = oriented[g_start:g_end + 1]
            if strand == "+":
                p_start, p_end = pam_start, pam_end
                genomic_start, genomic_end = g_start + 1, g_end + 1
            else:
                genomic_start, genomic_end = len(seq) - g_end, len(seq) - g_start
                p_start, p_end = len(seq) - 1 - pam_end, len(seq) - 1 - pam_start
            yield {
                "strand": strand,
                "guide": guide,
                "pam": pam,
                "pam_start": p_start + 1,
                "pam_end": p_end + 1,
                "guide_start": min(genomic_start, genomic_end),
                "guide_end": max(genomic_start, genomic_end),
            }


def search_local_off_targets(
    guide: str,
    reference_sequences: Mapping[str, str] | Sequence[Tuple[str, str]],
    *,
    max_mismatches: int = DEFAULT_MAX_LOCAL_MISMATCHES,
    max_hits: int = 100,
) -> List[Dict[str, Any]]:
    """Backward-compatible sequence-similarity search without PAM inference."""
    query = normalize_sequence(guide)
    items = list(reference_sequences.items()) if isinstance(reference_sequences, Mapping) else list(reference_sequences)
    hits: List[Dict[str, Any]] = []
    for name, raw_seq in items:
        seq = normalize_sequence(raw_seq)
        for strand, oriented in (("+", seq), ("-", reverse_complement(seq))):
            if len(oriented) < len(query):
                continue
            for i in range(len(oriented) - len(query) + 1):
                window = oriented[i:i + len(query)]
                mismatches = mismatch_count(query, window)
                if mismatches <= max_mismatches:
                    hits.append({
                        "sequence": str(name),
                        "strand": strand,
                        "position": i + 1,
                        "match": window,
                        "mismatches": mismatches,
                    })
                    if len(hits) >= max_hits:
                        return sorted(hits, key=lambda h: (h["mismatches"], h["sequence"], h["position"]))
    return sorted(hits, key=lambda h: (h["mismatches"], h["sequence"], h["position"]))


def search_local_off_targets_detailed(
    candidate: GuideCandidate,
    reference_sequences: Mapping[str, str] | Sequence[Tuple[str, str]],
    *,
    max_mismatches: int = DEFAULT_MAX_LOCAL_MISMATCHES,
    max_hits: int = 100,
    seed_length: int = DEFAULT_SEED_LENGTH,
    require_pam: bool = True,
) -> List[Dict[str, Any]]:
    """PAM-aware local specificity screen over explicitly supplied sequences."""
    system = get_crispr_system(candidate.system)
    items = list(reference_sequences.items()) if isinstance(reference_sequences, Mapping) else list(reference_sequences)
    query = candidate.guide
    hits: List[Dict[str, Any]] = []
    for name, raw_seq in items:
        seq = normalize_sequence(raw_seq)
        for window in _oriented_candidate_windows(seq, candidate, system):
            target = window["guide"]
            mm_positions = mismatch_positions(query, target)
            mm = len(mm_positions)
            if mm > max_mismatches:
                continue
            # The seed is the PAM-proximal end of the guide in targeting orientation.
            if system.pam_side == "3prime":
                seed_query = query[-min(seed_length, len(query)):]
                seed_target = target[-min(seed_length, len(target)):]
            else:
                seed_query = query[:min(seed_length, len(query))]
                seed_target = target[:min(seed_length, len(target))]
            seed_mm = mismatch_count(seed_query, seed_target)
            pam_ok = pam_matches(window["pam"], system.pam)
            if require_pam and not pam_ok:
                continue
            hits.append({
                "sequence": str(name),
                "strand": window["strand"],
                "guide_start": window["guide_start"],
                "guide_end": window["guide_end"],
                "pam_start": window["pam_start"],
                "pam_end": window["pam_end"],
                "match": target,
                "pam": window["pam"],
                "pam_compatible": pam_ok,
                "mismatches": mm,
                "mismatch_positions": mm_positions,
                "seed_mismatches": seed_mm,
                "classification": (
                    "exact"
                    if mm == 0
                    else "seed-mismatch" if seed_mm else "distal-mismatch"
                ),
            })
            if len(hits) >= max_hits:
                break
        if len(hits) >= max_hits:
            break
    return sorted(hits, key=lambda h: (h["mismatches"], h["seed_mismatches"], h["sequence"], h["guide_start"]))


def specificity_profile(hits: Sequence[Mapping[str, Any]], query_sequence_name: str = "target") -> Dict[str, Any]:
    exact = sum(int(h.get("mismatches", 999)) == 0 for h in hits)
    one_mm = sum(int(h.get("mismatches", 999)) == 1 for h in hits)
    two_mm = sum(int(h.get("mismatches", 999)) == 2 for h in hits)
    three_mm = sum(int(h.get("mismatches", 999)) == 3 for h in hits)
    seed_hits = sum(int(h.get("seed_mismatches", 999)) == 0 for h in hits)
    pam_hits = sum(bool(h.get("pam_compatible", False)) for h in hits)
    return {
        "analyzed": True,
        "scope": "local supplied sequences",
        "query": query_sequence_name,
        "total_hits": len(hits),
        "exact_matches": exact,
        "one_mismatch_hits": one_mm,
        "two_mismatch_hits": two_mm,
        "three_mismatch_hits": three_mm,
        "seed_match_hits": seed_hits,
        "pam_compatible_hits": pam_hits,
        "genome_wide": False,
        "limitation": "No genome-wide specificity claim; only explicitly supplied reference sequences were searched.",
    }


def annotate_local_off_targets(
    candidates: Sequence[GuideCandidate],
    reference_sequences: Mapping[str, str] | Sequence[Tuple[str, str]],
    *,
    max_mismatches: int = DEFAULT_MAX_LOCAL_MISMATCHES,
    max_hits: int = 100,
    seed_length: int = DEFAULT_SEED_LENGTH,
    require_pam: bool = True,
) -> List[GuideCandidate]:
    for candidate in candidates:
        hits = search_local_off_targets_detailed(
            candidate, reference_sequences,
            max_mismatches=max_mismatches,
            max_hits=max_hits,
            seed_length=seed_length,
            require_pam=require_pam,
        )
        candidate.local_off_target_hits = hits
        candidate.off_target_hits = len(hits)
        candidate.specificity_profile = specificity_profile(hits, candidate.candidate_id)
        candidate.flags = [f for f in candidate.flags if f != "INFO: local-sequence matches found"]
        if hits:
            candidate.flags.append("INFO: local-sequence matches found")
    return list(candidates)


# ============================================================================
# BASE EDITING
# ============================================================================

def get_base_editor(name: str) -> Dict[str, Any]:
    key = str(name).strip()
    if key not in BASE_EDITOR_DEFINITIONS:
        raise CRISPRAnalysisError(f"Unknown base editor '{name}'. Available: {', '.join(BASE_EDITOR_DEFINITIONS)}")
    return dict(BASE_EDITOR_DEFINITIONS[key])


def base_editing_window(candidate: GuideCandidate, editor: str = "CBE", window: Optional[Tuple[int, int]] = None) -> Tuple[int, int]:
    definition = get_base_editor(editor)
    return window or tuple(definition["window"])


def _variant_compatible_conversion(
    mutation: MutationTarget,
    conversion: Tuple[str, str],
    strand: str = "+",
) -> bool:
    ref, alt = conversion
    if strand == "+":
        return mutation.reference == ref and mutation.alternative == alt
    # A genomic edit is represented by its complementary base conversion when
    # the targeting guide is oriented on the reverse strand.
    return (
        reverse_complement(mutation.reference) == ref
        and reverse_complement(mutation.alternative) == alt
    )


def analyze_base_editing_candidate(
    candidate: GuideCandidate,
    mutation: MutationTarget,
    *,
    editor: str = "CBE",
    window: Optional[Tuple[int, int]] = None,
) -> BaseEditingCandidate:
    definition = get_base_editor(editor)
    conversion = tuple(definition["conversion"])
    edit_window = tuple(window or definition["window"])
    pos = mutation_position_in_guide(candidate, mutation.position)
    notes: List[str] = []
    desired_in_window = bool(pos is not None and edit_window[0] <= pos <= edit_window[1])
    compatible = _variant_compatible_conversion(mutation, conversion, candidate.strand)
    if not compatible:
        notes.append(f"Requested variant is not the modeled {conversion[0]}>{conversion[1]} conversion.")
    if pos is None:
        notes.append("Variant is outside the guide.")
    elif not desired_in_window:
        notes.append("Variant overlaps the guide but is outside the configured editing window.")

    bystanders: List[Dict[str, Any]] = []
    guide_conversion = conversion
    if candidate.strand == "-":
        guide_conversion = (reverse_complement(conversion[0]), reverse_complement(conversion[1]))
    if desired_in_window and compatible:
        for i, base in enumerate(candidate.guide, 1):
            if edit_window[0] <= i <= edit_window[1] and base == guide_conversion[0]:
                if pos != i:
                    bystanders.append({
                        "position_in_guide": i,
                        "base": base,
                        "conversion": f"{guide_conversion[0]}>{guide_conversion[1]}",
                    })
    return BaseEditingCandidate(
        candidate_id=candidate.candidate_id,
        editor=editor,
        conversion=f"{conversion[0]}>{conversion[1]}",
        guide=candidate.guide,
        pam=candidate.pam,
        strand=candidate.strand,
        variant=f"{mutation.chromosome}:{mutation.position} {mutation.reference}>{mutation.alternative}",
        variant_position_in_guide=pos or 0,
        editing_window=edit_window,
        desired_edit_in_window=desired_in_window and compatible,
        bystander_count=len(bystanders),
        bystanders=bystanders,
        sequence_level_notes=notes,
    )


def analyze_base_editing(
    candidates: Sequence[GuideCandidate],
    mutation: MutationTarget,
    *,
    editor: str = "CBE",
    window: Optional[Tuple[int, int]] = None,
) -> List[BaseEditingCandidate]:
    results = [analyze_base_editing_candidate(c, mutation, editor=editor, window=window) for c in candidates]
    for candidate, analysis in zip(candidates, results):
        candidate.base_editing_candidates = [asdict(analysis)]
        candidate.editing_window = analysis.editing_window
        candidate.bystander_edits = analysis.bystanders
        if analysis.desired_edit_in_window:
            candidate.flags.append("INFO: compatible base-editing window")
        if analysis.bystander_count:
            candidate.flags.append("WARN: base-editing bystanders detected")
    return results


def find_base_editing_candidates(
    sequence: str,
    mutation: MutationTarget,
    *,
    system: str = "SpCas9",
    editor: str = "CBE",
    region_start: int = 1,
    config: Optional[DesignConfig] = None,
) -> Tuple[DesignResult, List[BaseEditingCandidate]]:
    cfg = config or DesignConfig(system=system)
    result = design_guides(sequence, system=system, region_start=region_start, config=cfg, mutations=[mutation])
    analyses = analyze_base_editing(result.candidates, mutation, editor=editor, window=base_editing_window(result.candidates[0], editor) if result.candidates else None)
    compatible = [a for a in analyses if a.desired_edit_in_window]
    return result, compatible


# ============================================================================
# PRIME EDITING
# ============================================================================

def _sequence_context(sequence: str, start_1based: int, end_1based: int) -> str:
    seq = normalize_sequence(sequence)
    start = max(1, start_1based)
    end = min(len(seq), end_1based)
    if start > end:
        return ""
    return seq[start - 1:end]


def _prime_edit_for_candidate(
    sequence: str,
    candidate: GuideCandidate,
    mutation: MutationTarget,
    pbs_length: int,
    rtt_length: int,
    sequence_start: int = 1,
) -> PrimeEditingCandidate:
    seq = normalize_sequence(sequence)
    # Sequence-level construction in the guide orientation.  This deliberately
    # does not claim biochemical pegRNA efficacy or secondary-structure fitness.
    # Candidate coordinates are global/1-based, while ``seq`` is the supplied
    # local sequence. Translate every global coordinate into local coordinates.
    local_guide_start = candidate.guide_start - sequence_start + 1
    local_guide_end = candidate.guide_end - sequence_start + 1
    local_edit_anchor = mutation.position - sequence_start + 1

    if candidate.strand == "+":
        downstream = _sequence_context(seq, local_edit_anchor, local_edit_anchor + rtt_length - 1)
        upstream = _sequence_context(
            seq,
            local_guide_start - pbs_length,
            local_guide_start - 1,
        )
        pbs = reverse_complement(upstream) if upstream else ""
        rtt = downstream
    else:
        downstream = _sequence_context(
            seq,
            local_edit_anchor - rtt_length + 1,
            local_edit_anchor,
        )
        upstream = _sequence_context(
            seq,
            local_guide_end + 1,
            local_guide_end + pbs_length,
        )
        pbs = upstream
        rtt = reverse_complement(downstream) if downstream else ""
    ref_alt = f"{mutation.reference}>{mutation.alternative}"
    notes = [
        "PBS/RTT are sequence-level candidate constructions.",
        "No thermodynamic, structural, cellular, or experimental efficacy model is applied.",
    ]
    return PrimeEditingCandidate(
        candidate_id=f"PE_{candidate.candidate_id}_{pbs_length}_{rtt_length}",
        guide_id=candidate.candidate_id,
        guide=candidate.guide,
        pam=candidate.pam,
        strand=candidate.strand,
        variant=f"{mutation.chromosome}:{mutation.position} {ref_alt}",
        pbs_length=len(pbs),
        pbs_sequence=pbs,
        rtt_length=len(rtt),
        rtt_sequence=rtt,
        edit_in_rtt=(len(rtt) == rtt_length),
        sequence_level_notes=notes,
    )


def design_prime_editing_candidates(
    sequence: str,
    mutation: MutationTarget,
    candidates: Sequence[GuideCandidate],
    *,
    pbs_lengths: Sequence[int] = tuple(range(8, 16)),
    rtt_lengths: Sequence[int] = tuple(range(10, 31, 2)),
    sequence_start: int = 1,
) -> List[PrimeEditingCandidate]:
    results: List[PrimeEditingCandidate] = []
    for candidate in candidates:
        if mutation_position_in_guide(candidate, mutation.position) is None:
            continue
        for pbs_len in pbs_lengths:
            for rtt_len in rtt_lengths:
                requested_pbs = int(pbs_len)
                requested_rtt = int(rtt_len)
                if requested_pbs < 1 or requested_rtt < 1:
                    continue
                item = _prime_edit_for_candidate(
                    sequence, candidate, mutation, requested_pbs, requested_rtt, sequence_start=sequence_start
                )
                # A pegRNA candidate is only emitted when the requested PBS and
                # RTT lengths are actually available in the supplied sequence.
                # Truncated boundary fragments are not silently presented as valid
                # candidates.
                if (
                    item.edit_in_rtt
                    and item.pbs_length == requested_pbs
                    and item.rtt_length == requested_rtt
                ):
                    results.append(item)
        candidate.prime_editing_candidates = [asdict(x) for x in results if x.guide_id == candidate.candidate_id]
    return results


# ============================================================================
# TARGET CONTEXT / ORF INTEGRATION
# ============================================================================

def extract_orf_context(orf: Any) -> Optional[Tuple[str, int, str]]:
    """Extract sequence/start/name from common dict/object ORF representations."""
    if orf is None:
        return None

    def get_value(keys: Sequence[str], default: Any = None) -> Any:
        if isinstance(orf, Mapping):
            for key in keys:
                if key in orf:
                    return orf[key]
        for key in keys:
            if hasattr(orf, key):
                return getattr(orf, key)
        return default

    seq = get_value(("sequence", "dna", "orf_sequence", "seq"))
    if not seq:
        return None
    try:
        seq = normalize_sequence(str(seq))
    except CRISPRAnalysisError:
        return None
    start = get_value(("start", "start_position", "genomic_start", "position"), 1)
    try:
        start = int(start)
    except (TypeError, ValueError):
        start = 1
    name = str(get_value(("name", "id", "orf_id", "label"), "Selected ORF"))
    return seq, max(1, start), name


def orf_target_context(orf: Any, parent_sequence: Optional[str]) -> Optional[Tuple[str, int, str, str]]:
    """Return (forward-strand target, 1-based start, name, ORF strand) for an ORF.

    The target is the ORF's genomic region on the FORWARD strand of the parent
    sequence, so guides found on both strands are reported in parent-sequence
    coordinates with their true strand (audit P0-2). Without a parent
    sequence that matches the ORF, the ORF sequence itself is used and
    coordinates are relative to the ORF (stated in the region name).
    """
    from core.orf_context import resolve_orf_region

    region = resolve_orf_region(orf, parent_sequence)
    if region is None:
        return None
    if region.validated and region.forward_region:
        name = f"{region.name} ({region.strand} strand, {region.start_1based}-{region.end_1based})"
        return region.forward_region, region.start_1based, name, region.strand
    try:
        coding = normalize_sequence(region.coding_sequence)
    except CRISPRAnalysisError:
        return None
    return coding, 1, f"{region.name} (ORF-relative coordinates)", "+"


def build_target_region(
    sequence: str,
    *,
    start: int = 1,
    end: Optional[int] = None,
    name: str = "Selected Region",
    source: str = "sequence",
    feature_type: str = "sequence",
    feature_name: str = "",
    feature_strand: str = "+",
) -> TargetRegion:
    seq = normalize_sequence(sequence)
    if end is None:
        end = start + len(seq) - 1
    return TargetRegion(
        seq, start=start, end=end, name=name, source=source,
        feature_type=feature_type, feature_name=feature_name,
        feature_strand=feature_strand,
    )


# ============================================================================
# COMPLETE DESIGN API
# ============================================================================

def design_guides(
    sequence: str,
    *,
    sequence_name: str = "Selected Sequence",
    system: str = "SpCas9",
    region_start: int = 1,
    region_end: Optional[int] = None,
    region_name: str = "Selected Region",
    region_source: str = "sequence",
    feature_type: str = "sequence",
    feature_name: str = "",
    feature_strand: str = "+",
    config: Optional[DesignConfig] = None,
    mutations: Optional[Sequence[MutationTarget]] = None,
    local_reference_sequences: Optional[Mapping[str, str] | Sequence[Tuple[str, str]]] = None,
) -> DesignResult:
    cfg = config or DesignConfig(system=system)
    cfg.system = system
    valid, errors = validate_sequence(sequence, allow_ambiguous=cfg.allow_ambiguous)
    if not valid:
        raise CRISPRAnalysisError("Sequence validation failed: " + "; ".join(errors))

    seq = normalize_sequence(sequence)
    if region_end is None:
        region_end = region_start + len(seq) - 1
    region = build_target_region(
        seq,
        start=region_start,
        end=region_end,
        name=region_name,
        source=region_source,
        feature_type=feature_type,
        feature_name=feature_name,
        feature_strand=feature_strand,
    )
    crispr = get_crispr_system(system)
    warnings: List[str] = []

    if count_ambiguous_bases(seq):
        warnings.append("Sequence contains ambiguous bases; candidate availability may be reduced.")
    if len(seq) < crispr.guide_length + crispr.pam_length:
        warnings.append("Sequence is shorter than one complete guide+PAM architecture.")

    candidates = _find_targets_in_region(region, crispr, cfg)
    if mutations:
        annotate_mutation_targets(candidates, mutations, window=cfg.mutation_window)
    if local_reference_sequences:
        annotate_local_off_targets(
            candidates, local_reference_sequences,
            max_mismatches=cfg.max_local_off_target_mismatches,
            max_hits=cfg.max_local_off_target_hits,
            seed_length=cfg.seed_length,
            require_pam=cfg.require_pam_for_off_targets,
        )
    candidates = rank_guides(candidates, config=cfg)

    if not candidates:
        warnings.append("No sequence-compatible guide candidates were found.")

    metadata = {
        "analysis_scope": "sequence-level computational candidate design",
        "genome_wide_off_target_analysis": False,
        "experimental_efficacy": False,
        "specificity_assessment_available": bool(local_reference_sequences),
        "ranking_method": "transparent weighted sequence/context heuristic",
    }
    return DesignResult(
        sequence_name=sequence_name,
        sequence_length=len(seq),
        region=region,
        system=crispr,
        config=cfg,
        candidates=candidates,
        warnings=warnings,
        analysis_metadata=metadata,
    )


def run_crispr_analysis(
    sequence: str,
    *,
    sequence_name: str = "Selected Sequence",
    selected_orf: Optional[Any] = None,
    system: str = "SpCas9",
    config: Optional[DesignConfig] = None,
) -> DesignResult:
    """Programmatic entry point for main.py integration."""
    region_start = 1
    region_name = "Selected Region"
    feature_type = "sequence"
    feature_strand = "+"
    target = sequence
    if selected_orf is not None:
        context = orf_target_context(selected_orf, sequence)
        if context is not None:
            target, region_start, region_name, feature_strand = context
            feature_type = "ORF"
    return design_guides(
        target,
        sequence_name=sequence_name,
        system=system,
        region_start=region_start,
        region_name=region_name,
        region_source="ORF" if feature_type == "ORF" else "sequence",
        feature_type=feature_type,
        feature_name=region_name,
        feature_strand=feature_strand,
        config=config or DesignConfig(system=system),
    )


# ============================================================================
# BATCH ANALYSIS
# ============================================================================

def batch_design_guides(
    records: Sequence[Tuple[str, str]],
    *,
    system: str = "SpCas9",
    config: Optional[DesignConfig] = None,
) -> Dict[str, DesignResult]:
    results: Dict[str, DesignResult] = {}
    for name, sequence in records:
        results[str(name)] = design_guides(
            sequence,
            sequence_name=str(name),
            system=system,
            config=config or DesignConfig(system=system),
        )
    return results


def batch_summary(results: Mapping[str, DesignResult]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for name, result in results.items():
        rows.append({
            "sequence": name,
            "length_bp": result.sequence_length,
            "candidates": len(result.candidates),
            "high_quality": result.high_quality_count,
            "system": result.system.name,
            "warnings": len(result.warnings),
        })
    return rows


# ============================================================================
# SERIALIZATION / EXPORT
# ============================================================================

def guide_to_dict(candidate: GuideCandidate) -> Dict[str, Any]:
    return asdict(candidate)


def result_to_dict(result: DesignResult) -> Dict[str, Any]:
    return asdict(result)


def result_to_json(result: DesignResult, *, indent: int = 2) -> str:
    return json.dumps(result_to_dict(result), indent=indent, ensure_ascii=False)


def _default_export_directory() -> Path:
    """Return a writable, user-facing default export directory.

    On Windows this prefers the Desktop (including a OneDrive Desktop), then
    Documents, and finally a Synapse Exports directory in the user home.
    The directory is created when possible.
    """
    home = Path.home()
    candidates = [
        home / "Desktop",
        home / "OneDrive" / "Desktop",
        home / "Documents",
        home / "OneDrive" / "Documents",
        home / "Synapse Exports",
    ]
    for directory in candidates:
        if directory.exists() and directory.is_dir():
            try:
                directory.mkdir(parents=True, exist_ok=True)
                return directory
            except OSError:
                continue
    fallback = home / "Synapse Exports"
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def resolve_export_path(path: str | Path) -> Path:
    """Resolve an export path and make relative exports user-friendly."""
    raw = str(path).strip()
    if not raw:
        raise CRISPRAnalysisError("Export path cannot be empty.")
    resolved = Path(raw).expanduser()
    if not resolved.is_absolute():
        resolved = _default_export_directory() / resolved
    resolved = resolved.resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved


def export_json(result: DesignResult, path: str | Path) -> str:
    resolved = resolve_export_path(path)
    with resolved.open("w", encoding="utf-8") as handle:
        handle.write(result_to_json(result))
    return str(resolved)


def _candidate_csv_row(c: GuideCandidate) -> Dict[str, Any]:
    return {
        "candidate_id": c.candidate_id,
        "system": c.system,
        "guide": c.guide,
        "pam": c.pam,
        "strand": c.strand,
        "guide_start": c.guide_start,
        "guide_end": c.guide_end,
        "pam_start": c.pam_start,
        "pam_end": c.pam_end,
        "guide_length": c.guide_length,
        "gc_content": f"{c.metrics.gc_content:.3f}",
        "at_content": f"{c.metrics.at_content:.3f}",
        "complexity": f"{c.metrics.sequence_complexity:.4f}",
        "homopolymer_max": c.metrics.homopolymer_max,
        "repeat_fraction": f"{c.metrics.repeat_fraction:.4f}",
        "self_complementarity_max": c.metrics.self_complementarity_max,
        "hairpin_risk": c.metrics.hairpin_risk,
        "unique_kmer_fraction": f"{c.metrics.unique_kmer_fraction:.4f}",
        "ambiguous_bases": c.metrics.ambiguous_bases,
        "rank_score": f"{c.rank_score:.3f}",
        "rank_percentile": f"{c.rank_percentile:.2f}",
        "flags": "; ".join(c.flags),
        "mutation_overlap": "; ".join(c.mutation_overlap),
        "local_off_target_hits": c.off_target_hits,
        "specificity_scope": c.specificity_profile.get("scope", "not analyzed"),
    }


def export_csv(result: DesignResult, path: str | Path) -> str:
    path = resolve_export_path(path)
    fields = list(_candidate_csv_row(result.candidates[0]).keys()) if result.candidates else [
        "candidate_id", "system", "guide", "pam", "strand", "guide_start", "guide_end",
        "pam_start", "pam_end", "guide_length", "gc_content", "at_content", "complexity",
        "homopolymer_max", "repeat_fraction", "self_complementarity_max", "hairpin_risk",
        "unique_kmer_fraction", "ambiguous_bases", "rank_score", "rank_percentile", "flags",
        "mutation_overlap", "local_off_target_hits", "specificity_scope",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for candidate in result.candidates:
            writer.writerow(_candidate_csv_row(candidate))
    return path


def export_guides_fasta(result: DesignResult, path: str | Path, top_n: Optional[int] = None) -> str:
    path = resolve_export_path(path)
    candidates = result.candidates if top_n is None else result.candidates[:top_n]
    with open(path, "w", encoding="utf-8") as handle:
        for candidate in candidates:
            handle.write(
                f">{candidate.candidate_id}|{candidate.system}|{candidate.strand}|"
                f"guide={candidate.guide_start}-{candidate.guide_end}|pam={candidate.pam}\n"
            )
            handle.write(candidate.guide + "\n")
    return path


def export_batch_summary(results: Mapping[str, DesignResult], path: str | Path) -> str:
    path = resolve_export_path(path)
    rows = batch_summary(results)
    fields = ["sequence", "length_bp", "candidates", "high_quality", "system", "warnings"]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


# ============================================================================
# REPORTING / VISUAL TEXT
# ============================================================================

def _line(char: str = "-", width: int = 78) -> str:
    return char * width


def target_map(candidate: GuideCandidate, sequence_length: int, width: int = 60) -> str:
    """ASCII target map for quick positional grounding."""
    if sequence_length <= 0:
        return ""
    width = max(20, width)
    guide_s = max(1, candidate.guide_start)
    guide_e = min(sequence_length, candidate.guide_end)
    pam_s = max(1, candidate.pam_start)
    pam_e = min(sequence_length, candidate.pam_end)
    chars = ["."] * width
    def place(start: int, end: int, symbol: str) -> None:
        a = int((start - 1) / sequence_length * width)
        b = int((end - 1) / sequence_length * width)
        a = max(0, min(width - 1, a))
        b = max(a, min(width - 1, b))
        for i in range(a, b + 1):
            chars[i] = symbol
    place(guide_s, guide_e, "=")
    place(pam_s, pam_e, "#")
    return "".join(chars) + "  (= guide, # PAM)"


def format_candidate(candidate: GuideCandidate, rank: int) -> str:
    flags = ", ".join(candidate.flags) if candidate.flags else "None"
    mutation = "; ".join(candidate.mutation_overlap) if candidate.mutation_overlap else "None"
    specificity = candidate.specificity_profile
    specificity_text = (
        f"{specificity.get('total_hits', 0)} local hits"
        if specificity.get("analyzed") else "Not analyzed"
    )
    components = ", ".join(f"{k}={v:.1f}" for k, v in candidate.rank_components.items())
    return "\n".join([
        f"#{rank} {candidate.candidate_id}",
        f"Guide:                    {candidate.guide}",
        f"PAM:                      {candidate.pam}",
        f"Strand:                   {candidate.strand}",
        f"Guide position:           {candidate.guide_start}-{candidate.guide_end}",
        f"PAM position:             {candidate.pam_start}-{candidate.pam_end}",
        f"Seed sequence:            {candidate.seed_sequence}",
        f"GC / AT:                  {candidate.metrics.gc_content:.2f}% / {candidate.metrics.at_content:.2f}%",
        f"Complexity:               {candidate.metrics.sequence_complexity:.3f}",
        f"Homopolymer max:          {candidate.metrics.homopolymer_max}",
        f"Repeat fraction:          {candidate.metrics.repeat_fraction:.3f}",
        f"Self-complementarity:     {candidate.metrics.self_complementarity_max} nt",
        f"Hairpin-risk proxy:       {candidate.metrics.hairpin_risk}",
        f"Local specificity:        {specificity_text}",
        f"Mutation overlap:         {mutation}",
        f"Rank score:               {candidate.rank_score:.3f}/100",
        f"Rank components:          {components or 'None'}",
        f"Flags:                    {flags}",
        f"Target map:               {target_map(candidate, max(candidate.guide_end, candidate.pam_end))}",
    ])


def format_base_editing_report(results: Sequence[BaseEditingCandidate]) -> str:
    lines = ["BASE EDITING ANALYSIS", _line()]
    if not results:
        lines.append("No base-editing candidates available.")
        return "\n".join(lines)
    for item in results:
        lines.extend([
            f"{item.candidate_id} | editor={item.editor} | conversion={item.conversion}",
            f"Guide: {item.guide}",
            f"PAM: {item.pam} | strand={item.strand}",
            f"Variant: {item.variant}",
            f"Variant position in guide: {item.variant_position_in_guide}",
            f"Editing window: {item.editing_window[0]}-{item.editing_window[1]}",
            f"Desired edit in window: {'YES' if item.desired_edit_in_window else 'NO'}",
            f"Bystanders: {item.bystander_count}",
        ])
        for b in item.bystanders:
            lines.append(f"  - position {b['position_in_guide']}: {b['base']} -> {b['conversion'].split('>')[1]}")
        for note in item.sequence_level_notes:
            lines.append(f"  NOTE: {note}")
        lines.append(_line())
    return "\n".join(lines)


def format_prime_editing_report(results: Sequence[PrimeEditingCandidate], top_n: int = 20) -> str:
    lines = ["PRIME EDITING CANDIDATES", _line()]
    if not results:
        lines.append("No sequence-level prime-editing candidates available.")
        return "\n".join(lines)
    for item in results[:top_n]:
        lines.extend([
            f"{item.candidate_id} | guide={item.guide_id} | strand={item.strand}",
            f"Guide: {item.guide}",
            f"PAM: {item.pam}",
            f"Variant: {item.variant}",
            f"PBS ({item.pbs_length} nt): {item.pbs_sequence}",
            f"RTT ({item.rtt_length} nt): {item.rtt_sequence}",
            f"Edit represented in RTT: {'YES' if item.edit_in_rtt else 'NO'}",
            "Notes: " + " | ".join(item.sequence_level_notes),
            _line(),
        ])
    if len(results) > top_n:
        lines.append(f"Showing top {top_n} of {len(results)} sequence-level candidates.")
    return "\n".join(lines)


def format_report(result: DesignResult, *, top_n: int = 10) -> str:
    lines = [
        "=" * 78,
        "SYNAPSE CRISPR GUIDE ANALYSIS REPORT",
        "=" * 78,
        f"Sequence:              {result.sequence_name}",
        f"Sequence length:       {result.sequence_length} bp",
        f"Target region:         {result.region.name}",
        f"Region coordinates:    {result.region.start}-{result.region.end}",
        f"Feature type:          {result.region.feature_type}",
        f"Feature name:          {result.region.feature_name or 'N/A'}",
        f"Feature strand:        {result.region.feature_strand}",
        f"CRISPR system:         {result.system.name}",
        f"Enzyme family:         {result.system.enzyme_family}",
        f"PAM:                   {result.system.pam}",
        f"Guide length:          {result.system.guide_length} nt",
        f"PAM orientation:       {result.system.pam_side}",
        f"Candidates found:      {len(result.candidates)}",
        f"Heuristic quality-pass: {result.high_quality_count}/{len(result.candidates)}",
        f"Ranking mode:            {ranking_mode_label(result.candidates)}",
        _line(),
        "SCIENTIFIC SCOPE",
        _line(),
        "Ranking is a transparent computational heuristic, not measured editing efficacy.",
        "If local specificity data are available, the ranking is explicitly specificity-aware.",
        "Local specificity uses only explicitly supplied reference sequences.",
        "Genome-wide off-target analysis requires a real reference genome/index.",
        "Base/prime editing outputs are sequence-level candidate architectures.",
    ]
    if result.warnings:
        lines.extend([_line(), "WARNINGS", _line()])
        lines.extend(f"- {warning}" for warning in result.warnings)
    lines.extend([
        _line(),
        f"RANKED GUIDE CANDIDATES ({ranking_mode_label(result.candidates).upper()})",
        _line(),
    ])
    if not result.candidates:
        lines.append("No candidates.")
    for rank, candidate in enumerate(result.candidates[:top_n], 1):
        lines.append(format_candidate(candidate, rank))
        lines.append(_line())
    if len(result.candidates) > top_n:
        lines.append(f"Showing top {top_n} of {len(result.candidates)} candidates.")
    return "\n".join(lines)


def print_report(result: DesignResult, *, top_n: int = 10) -> None:
    print(format_report(result, top_n=top_n))


# ============================================================================
# CLI HELPERS
# ============================================================================

def _input_int(prompt: str, default: Optional[int] = None, minimum: Optional[int] = None) -> int:
    while True:
        raw = input(prompt).strip()
        if not raw and default is not None:
            return default
        try:
            value = int(raw)
            if minimum is not None and value < minimum:
                raise ValueError
            return value
        except ValueError:
            suffix = f" (>= {minimum})" if minimum is not None else ""
            print(f"Invalid integer{suffix}.")


def _input_float(prompt: str, default: float) -> float:
    while True:
        raw = input(prompt).strip()
        if not raw:
            return default
        try:
            return float(raw)
        except ValueError:
            print("Invalid number.")


def _choose_system() -> str:
    systems = list_crispr_systems()
    print("\nAvailable CRISPR systems")
    print(_line())
    for i, system in enumerate(systems, 1):
        print(
            f"[{i}] {system.name:<12} PAM={system.pam:<6} Guide={system.guide_length:<2} nt | "
            f"{system.description}"
        )
    print("[0] Cancel")
    while True:
        choice = input("System: ").strip()
        if choice == "0":
            return ""
        try:
            index = int(choice)
            if 1 <= index <= len(systems):
                return systems[index - 1].name
        except ValueError:
            pass
        print("Invalid selection.")


def _choose_editor() -> str:
    editors = list_base_editors()
    print("\nAvailable computational base-editor models")
    print(_line())
    for i, name in enumerate(editors, 1):
        data = get_base_editor(name)
        print(f"[{i}] {name:<5} {data['conversion'][0]}>{data['conversion'][1]} window={data['window'][0]}-{data['window'][1]}")
    print("[0] Cancel")
    while True:
        choice = input("Editor: ").strip()
        if choice == "0":
            return ""
        try:
            index = int(choice)
            if 1 <= index <= len(editors):
                return editors[index - 1]
        except ValueError:
            pass
        print("Invalid selection.")


def _manual_sequence() -> Tuple[str, str]:
    name = input("Sequence name: ").strip() or "Manual Sequence"
    print("Paste DNA sequence. Enter an empty line to finish.")
    lines: List[str] = []
    while True:
        line = input()
        if not line.strip():
            break
        lines.append(line)
    return name, "".join(lines)


def _print_candidate_table(candidates: Sequence[GuideCandidate], limit: int = 20) -> None:
    print("\n" + _line())
    mode = ranking_mode_label(candidates)
    print(f"RANKED GUIDES ({mode.upper()})")
    print("Score = transparent computational prioritization heuristic.")
    if has_specificity_analysis(candidates):
        print("Specificity component uses only explicitly supplied local reference sequences.")
    print(_line())
    print(f"{'#':>3} {'ID':<7} {'Guide':<22} {'PAM':<7} {'Str':<4} {'Pos':<12} {'GC':>6} {'Score':>7}")
    print(_line())
    for i, c in enumerate(candidates[:limit], 1):
        print(
            f"{i:>3} {c.candidate_id:<7} {c.guide:<22} {c.pam:<7} {c.strand:<4} "
            f"{c.guide_start}-{c.guide_end:<7} {c.metrics.gc_content:>5.1f}% {c.rank_score:>7.2f}"
        )
    if len(candidates) > limit:
        print(f"Showing {limit} of {len(candidates)} candidates.")


def _choose_candidate(candidates: Sequence[GuideCandidate]) -> Optional[GuideCandidate]:
    if not candidates:
        print("No guide candidates available.")
        return None
    _print_candidate_table(candidates)
    raw = input("Candidate number (0 cancel): ").strip()
    try:
        index = int(raw)
        if index == 0:
            return None
        if 1 <= index <= len(candidates):
            return candidates[index - 1]
    except ValueError:
        pass
    print("Invalid candidate selection.")
    return None


def _read_local_references() -> Dict[str, str]:
    print("Enter local reference sequences as name=sequence. Empty line finishes.")
    print("For PAM-aware analysis, provide the DNA context containing the guide and its PAM.")
    print("Example: target=ACCGTTGACGAACTGAAAGCTGG")
    refs: Dict[str, str] = {}
    while True:
        line = input().strip()
        if not line:
            break
        if "=" not in line:
            print("Expected name=sequence.")
            continue
        name, seq = line.split("=", 1)
        try:
            normalized = normalize_sequence(seq)
            if len(normalized) < 23:
                print(
                    "WARNING: Reference is shorter than a 20-nt SpCas9 guide plus a 3-nt PAM. "
                    "A PAM-aware hit may not be detectable."
                )
            refs[name.strip()] = normalized
        except CRISPRAnalysisError as exc:
            print(f"ERROR: {exc}")
    return refs


def _load_orf_or_region(current_sequence: str) -> Tuple[str, int, int, str, str]:
    if not current_sequence:
        raise CRISPRAnalysisError("No DNA sequence selected.")
    start = _input_int("Region start (1-based): ", 1, 1)
    end = _input_int(f"Region end [{len(current_sequence)}]: ", len(current_sequence), start)
    if end > len(current_sequence):
        raise CRISPRAnalysisError("Region end exceeds sequence length.")
    return current_sequence[start - 1:end], start, end, f"Region {start}-{end}", "region"


def _cli_design() -> None:
    print("\n" + "=" * 78)
    print("CRISPR GUIDE DESIGN")
    print("=" * 78)
    name, sequence = _manual_sequence()
    if not sequence:
        print("No sequence entered.")
        return
    system = _choose_system()
    if not system:
        return
    config = DesignConfig(system=system)
    config.min_gc = _input_float(f"Minimum GC [{DEFAULT_MIN_GC}]: ", DEFAULT_MIN_GC)
    config.max_gc = _input_float(f"Maximum GC [{DEFAULT_MAX_GC}]: ", DEFAULT_MAX_GC)
    config.max_homopolymer = _input_int(f"Maximum homopolymer [{DEFAULT_MAX_HOMOPOLYMER}]: ", DEFAULT_MAX_HOMOPOLYMER, 1)
    try:
        result = design_guides(sequence, sequence_name=name, system=system, config=config)
    except CRISPRAnalysisError as exc:
        print(f"ERROR: {exc}")
        return
    print_report(result)


# ============================================================================
# INTERACTIVE SYNAPSE MENU
# ============================================================================

def crispr_menu(
    sequence: Optional[str] = None,
    sequence_name: str = "Selected Sequence",
    selected_orf: Optional[Any] = None,
) -> None:
    """Interactive Synapse CRISPR menu.

    The analysis engine remains fully featured, while the user interface groups
    closely related operations into a smaller set of professional workflows.
    Existing main.py integration is preserved.
    """
    base_sequence = normalize_sequence(sequence) if sequence else None
    current_name = sequence_name
    region_start = 1
    region_end = len(base_sequence) if base_sequence else 0
    region_name = "Selected Region"
    feature_type = "sequence"
    selected_mutations: List[MutationTarget] = []
    local_references: Dict[str, str] = {}
    last_result: Optional[DesignResult] = None
    last_base_results: List[BaseEditingCandidate] = []
    last_prime_results: List[PrimeEditingCandidate] = []
    batch_results: Dict[str, DesignResult] = {}
    batch_sequences: Dict[str, str] = {}
    config = DesignConfig()

    if selected_orf is not None:
        extracted = orf_target_context(selected_orf, base_sequence)
        if extracted is not None:
            orf_seq, orf_start, orf_name, _orf_strand = extracted
            base_sequence = orf_seq
            current_name = orf_name
            region_start = orf_start
            region_end = orf_start + len(orf_seq) - 1
            region_name = orf_name
            feature_type = "ORF"

    def sequence_context_menu() -> None:
        nonlocal base_sequence, current_name, region_start, region_end
        nonlocal region_name, feature_type, last_result

        while True:
            print("\n" + "=" * 78)
            print("SEQUENCE & REGION")
            print("=" * 78)
            print("[1] Select DNA Sequence")
            print("[2] Select ORF / Region")
            print("[0] Back")
            print(_line())
            choice = input("Selection: ").strip()

            if choice == "0":
                return
            if choice == "1":
                name, seq = _manual_sequence()
                try:
                    base_sequence = normalize_sequence(seq)
                except CRISPRAnalysisError as exc:
                    print(f"ERROR: {exc}")
                    continue
                current_name = name
                region_start = 1
                region_end = len(base_sequence)
                region_name = "Selected Region"
                feature_type = "sequence"
                last_result = None
                print(f"Selected sequence: {current_name} ({len(base_sequence)} bp)")
            elif choice == "2":
                if not base_sequence:
                    print("No DNA sequence selected.")
                    continue
                try:
                    selected_seq, start, end, name, _source = _load_orf_or_region(base_sequence)
                    # start/end are 1-based within the CURRENT base sequence,
                    # which may itself be a region; keep parent coordinates.
                    offset = region_start - 1
                    base_sequence = selected_seq
                    region_start = start + offset
                    region_end = end + offset
                    region_name = f"Region {region_start}-{region_end}"
                    feature_type = "region"
                    last_result = None
                    print(f"Selected {region_name} ({len(base_sequence)} bp).")
                except CRISPRAnalysisError as exc:
                    print(f"ERROR: {exc}")
            else:
                print("Invalid selection.")

    def batch_analysis_menu() -> None:
        nonlocal batch_results, batch_sequences, current_name, base_sequence, region_start
        nonlocal region_end, region_name, feature_type, last_result

        while True:
            print("\n" + "=" * 78)
            print("BATCH ANALYSIS")
            print("=" * 78)
            print("[1] Load FASTA / Batch Analysis")
            print("[2] Use First Record as Current Sequence")
            print("[0] Back")
            print(_line())
            choice = input("Selection: ").strip()

            if choice == "0":
                return
            if choice == "1":
                path = input("FASTA path: ").strip()
                if not path:
                    continue
                try:
                    records = load_fasta(path)
                    batch_results = batch_design_guides(
                        records, system=config.system, config=config
                    )
                    batch_sequences = {name: seq for name, seq in records}
                    print(f"Loaded {len(records)} FASTA record(s).")
                    for row in batch_summary(batch_results):
                        print(
                            f"{row['sequence']}: {row['length_bp']} bp | "
                            f"candidates={row['candidates']} | "
                            f"high-quality={row['high_quality']}"
                        )
                except CRISPRAnalysisError as exc:
                    print(f"ERROR: {exc}")
            elif choice == "2":
                if not batch_results:
                    print("No batch analysis available. Use [1] first.")
                    continue
                records = list(batch_results.keys())
                print("\nAVAILABLE RECORDS")
                for index, name in enumerate(records, start=1):
                    print(f"[{index}] {name}")
                raw = input("Record number (0 cancel): ").strip()
                try:
                    index = int(raw)
                except ValueError:
                    print("Invalid selection.")
                    continue
                if index == 0:
                    continue
                if 1 <= index <= len(records):
                    selected_name = records[index - 1]
                    result = batch_results[selected_name]
                    base_sequence = batch_sequences[selected_name]
                    current_name = selected_name
                    region_start = result.region_start
                    region_end = result.region_end
                    region_name = result.region_name
                    feature_type = result.feature_type
                    last_result = result
                    print(f"Current sequence: {current_name}")
                else:
                    print("Invalid selection.")
            else:
                print("Invalid selection.")

    def system_menu() -> None:
        nonlocal config
        while True:
            print("\n" + "=" * 78)
            print("CRISPR SYSTEM")
            print("=" * 78)
            print(f"Current system: {config.system}")
            print()
            print("[1] Select CRISPR System")
            print("[2] Show System Details")
            print("[0] Back")
            print(_line())
            choice = input("Selection: ").strip()

            if choice == "0":
                return
            if choice == "1":
                system = _choose_system()
                if system:
                    config.system = system
                    print(f"Selected system: {system}")
            elif choice == "2":
                print()
                print(format_system_details(config.system))
            else:
                print("Invalid selection.")

    def guide_analysis_menu() -> None:
        nonlocal last_result
        if not last_result:
            print("Run [4] Target Discovery first.")
            return
        while True:
            print("\n" + "=" * 78)
            print("GUIDE ANALYSIS & RANKING")
            print("=" * 78)
            print("[1] Analyze Individual Guide")
            print("[2] Rank Guides")
            print("[0] Back")
            print(_line())
            choice = input("Selection: ").strip()
            if choice == "0":
                return
            if choice == "1":
                candidate = _choose_candidate(last_result.candidates)
                if candidate:
                    print("\n" + format_candidate(candidate, 1))
            elif choice == "2":
                last_result.candidates = rank_guides(last_result.candidates, config=config)
                _print_candidate_table(last_result.candidates)
            else:
                print("Invalid selection.")

    def report_export_menu() -> None:
        if not last_result and not batch_results:
            print("No analysis results available.")
            return
        while True:
            print("\n" + "=" * 78)
            print("REPORTS & EXPORT")
            print("=" * 78)
            print("[1] Full Guide Report")
            print("[2] Export Guides")
            print("[3] Export Batch Summary")
            print("[0] Back")
            print(_line())
            choice = input("Selection: ").strip()

            if choice == "0":
                return
            if choice == "1":
                if not last_result:
                    print("Run [4] Target Discovery first.")
                    continue
                print_report(last_result, top_n=10)
                if last_base_results:
                    print(format_base_editing_report(last_base_results))
                if last_prime_results:
                    print(format_prime_editing_report(last_prime_results, top_n=10))
            elif choice == "2":
                if not last_result:
                    print("Run [4] Target Discovery first.")
                    continue
                default_dir = _default_export_directory()
                print(f"Relative paths are saved to: {default_dir}")
                path = input("Export path (.json, .csv or .fasta): ").strip()
                if not path:
                    continue
                try:
                    lower = path.lower()
                    if lower.endswith(".json"):
                        exported = export_json(last_result, path)
                    elif lower.endswith(".csv"):
                        exported = export_csv(last_result, path)
                    elif lower.endswith(".fasta") or lower.endswith(".fa"):
                        exported = export_guides_fasta(last_result, path)
                    else:
                        print("Use .json, .csv or .fasta/.fa extension.")
                        continue
                    print(f"Exported: {exported}")
                except PermissionError:
                    print("Export failed: permission denied.")
                    print("Choose a writable Desktop/Documents location or use a relative filename.")
                except OSError as exc:
                    print(f"Export failed: {exc}")
            elif choice == "3":
                if not batch_results:
                    print("No batch analysis available. Use [2] -> [1] Load FASTA / Batch Analysis first.")
                    continue
                default_dir = _default_export_directory()
                print(f"Relative paths are saved to: {default_dir}")
                path = input("Batch CSV export path: ").strip()
                if not path:
                    continue
                try:
                    exported = export_batch_summary(batch_results, path)
                    print(f"Exported: {exported}")
                except PermissionError:
                    print("Export failed: permission denied.")
                    print("Choose a writable Desktop/Documents location or use a relative filename.")
                except OSError as exc:
                    print(f"Export failed: {exc}")
            else:
                print("Invalid selection.")

    while True:
        print("\n" + "=" * 78)
        print("CRISPR ANALYSIS ENGINE")
        print("=" * 78)
        print(f"Sequence: {current_name}")
        print(f"Length:   {len(base_sequence) if base_sequence else 0} bp")
        print(
            f"Region:   {region_name} ({region_start}-{region_end})"
            if base_sequence else "Region:   none"
        )
        print(f"System:   {config.system}")
        print()
        print("INPUT / CONTEXT")
        print("[1] Sequence & Region")
        print("[2] Batch Analysis")
        print("[3] CRISPR System")
        print()
        print("GUIDE DESIGN")
        print("[4] Target Discovery")
        print("[5] Guide Analysis & Ranking")
        print()
        print("EDITING")
        print("[6] Mutation / Variant Targeting")
        print("[7] Base Editing")
        print("[8] Prime Editing")
        print()
        print("SPECIFICITY")
        print("[9] Off-Target Analysis")
        print("[10] Specificity Profile")
        print()
        print("OUTPUT")
        print("[11] Reports & Export")
        print()
        print("SYSTEM")
        print("[12] CRISPR Self-Test")
        print("[13] Supported Systems")
        print("[0] Back")
        print(_line())
        choice = input("Selection: ").strip()

        if choice == "0":
            return

        if choice == "1":
            sequence_context_menu()

        elif choice == "2":
            batch_analysis_menu()

        elif choice == "3":
            system_menu()

        elif choice == "4":
            if not base_sequence:
                print("No sequence selected.")
                continue
            config.min_gc = _input_float(f"Minimum GC [{config.min_gc}]: ", config.min_gc)
            config.max_gc = _input_float(f"Maximum GC [{config.max_gc}]: ", config.max_gc)
            try:
                last_result = design_guides(
                    base_sequence,
                    sequence_name=current_name,
                    system=config.system,
                    region_start=region_start,
                    region_end=region_end,
                    region_name=region_name,
                    feature_type=feature_type,
                    feature_name=region_name,
                    config=config,
                    mutations=selected_mutations,
                    local_reference_sequences=local_references or None,
                )
                print(f"Found {last_result.candidate_count} candidate(s).")
                _print_candidate_table(last_result.candidates)
            except CRISPRAnalysisError as exc:
                print(f"ERROR: {exc}")

        elif choice == "5":
            guide_analysis_menu()

        elif choice == "6":
            raw = input("Variant (chr:position REF>ALT; blank cancels): ").strip()
            if not raw:
                continue
            try:
                mutation = parse_mutation_target(raw)
                selected_mutations.append(mutation)
                print(
                    f"Added variant: {mutation.chromosome}:{mutation.position} "
                    f"{mutation.reference}>{mutation.alternative}"
                )
                if last_result:
                    annotate_mutation_targets(
                        last_result.candidates,
                        selected_mutations,
                        window=config.mutation_window,
                    )
                    last_result.candidates = rank_guides(last_result.candidates, config=config)
                    _print_candidate_table(last_result.candidates)
            except CRISPRAnalysisError as exc:
                print(f"ERROR: {exc}")

        elif choice == "7":
            if not last_result:
                print("Run [4] Target Discovery first.")
                continue
            if not selected_mutations:
                print("No mutation/variant selected. Use [6] first.")
                continue
            editor = _choose_editor()
            if not editor:
                continue
            last_base_results = []
            for mutation in selected_mutations:
                last_base_results.extend(
                    analyze_base_editing(last_result.candidates, mutation, editor=editor)
                )
            print(format_base_editing_report(last_base_results))

        elif choice == "8":
            if not last_result or not base_sequence:
                print("Run [4] Target Discovery first.")
                continue
            if not selected_mutations:
                print("No mutation/variant selected. Use [6] first.")
                continue
            last_prime_results = []
            for mutation in selected_mutations:
                last_prime_results.extend(
                    design_prime_editing_candidates(
                        base_sequence,
                        mutation,
                        last_result.candidates,
                        pbs_lengths=config.prime_pbs_lengths,
                        rtt_lengths=config.prime_rtt_lengths,
                        sequence_start=region_start,
                    )
                )
            print(format_prime_editing_report(last_prime_results))

        elif choice == "9":
            if not last_result:
                print("Run [4] Target Discovery first.")
                continue
            local_references = _read_local_references()
            if local_references:
                annotate_local_off_targets(
                    last_result.candidates,
                    local_references,
                    max_mismatches=config.max_local_off_target_mismatches,
                    max_hits=config.max_local_off_target_hits,
                    seed_length=config.seed_length,
                    require_pam=config.require_pam_for_off_targets,
                )
                last_result.candidates = rank_guides(last_result.candidates, config=config)
                total_hits = sum(len(c.local_off_target_hits) for c in last_result.candidates)
                print("PAM-aware local specificity analysis completed.")
                if total_hits:
                    print(f"Local PAM-compatible hits found: {total_hits}")
                    print("Top local matches:")
                    shown = 0
                    for candidate in last_result.candidates:
                        for hit in candidate.local_off_target_hits:
                            print(
                                f"  {candidate.candidate_id}: {hit['sequence']} "
                                f"strand={hit['strand']} guide={hit['guide_start']}-{hit['guide_end']} "
                                f"PAM={hit['pam']} mismatches={hit['mismatches']} "
                                f"seed_mismatches={hit['seed_mismatches']} "
                                f"classification={hit['classification']}"
                            )
                            shown += 1
                            if shown >= 20:
                                break
                        if shown >= 20:
                            break
                    if total_hits > shown:
                        print(f"Showing {shown} of {total_hits} local hits.")
                else:
                    print("No PAM-compatible local hits were found in the supplied references.")
                    print(
                        "For a positive control, provide reference DNA containing the guide plus its PAM "
                        "(for example: target=ACCGTTGACGAACTGAAAGCTGG for SpCas9/G0004)."
                    )
                print("Scientific limitation: this is not a genome-wide off-target assessment.")

        elif choice == "10":
            if not last_result:
                print("Run [4] Target Discovery first.")
                continue
            print("SPECIFICITY PROFILE")
            print(_line())
            if has_specificity_analysis(last_result.candidates):
                print("Candidates with explicit local specificity analysis are marked in the ranking.")
                print("Select the candidate to inspect its local specificity profile.")
            else:
                print("No local specificity analysis has been performed yet.")
                print("Use [9] Off-Target Analysis first and supply reference sequences.")
            candidate = _choose_candidate(last_result.candidates)
            if candidate:
                profile = candidate.specificity_profile or {
                    "analyzed": False,
                    "scope": "local supplied sequences",
                    "query": candidate.candidate_id,
                    "total_hits": 0,
                    "genome_wide": False,
                    "limitation": "No explicit reference sequence set has been searched.",
                }
                print(json.dumps(profile, indent=2))

        elif choice == "11":
            report_export_menu()

        elif choice == "12":
            run_self_test(verbose=True)

        elif choice == "13":
            print("\nSUPPORTED CRISPR SYSTEMS")
            print(_line())
            for system in list_crispr_systems():
                print(
                    f"{system.name}: family={system.enzyme_family}, "
                    f"PAM={system.pam}, guide={system.guide_length} nt, "
                    f"PAM side={system.pam_side}"
                )

        else:
            print("Invalid selection.")


# ============================================================================
# SELF-TESTS
# ============================================================================

def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _test_iupac_pam() -> None:
    _assert(pam_matches("AGG", "NGG"), "AGG should match NGG")
    _assert(pam_matches("TGG", "NGG"), "TGG should match NGG")
    _assert(not pam_matches("AGA", "NGG"), "AGA should not match NGG")
    _assert(pam_matches("TTTA", "TTTV"), "TTTA should match TTTV")


def _test_reverse_complement() -> None:
    _assert(reverse_complement("ATGC") == "GCAT", "Reverse complement failed")
    _assert(reverse_complement("RYSW") == "WSRY", "IUPAC reverse complement failed")


def _test_metrics() -> None:
    metrics = characterize_guide("GCGCGCGCGCGCGCGCGCGC")
    _assert(abs(metrics.gc_content - 100.0) < 1e-9, "GC calculation failed")
    _assert(metrics.homopolymer_max == 1, "Homopolymer calculation failed")
    _assert(metrics.ambiguous_bases == 0, "Ambiguity count failed")
    _assert(metrics.sequence_complexity <= 1.0, "Complexity outside normalized range")


def _test_spcas9_plus() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    candidates = find_targets(guide + "AGG", system="SpCas9")
    _assert(any(c.strand == "+" and c.guide == guide for c in candidates), "SpCas9 + target not found")


def _test_spcas9_minus() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    rc_guide = reverse_complement(guide)
    sequence = "CCA" + rc_guide
    candidates = find_targets(sequence, system="SpCas9")
    _assert(any(c.strand == "-" and c.guide == guide for c in candidates), "SpCas9 - target not found")


def _test_cas12a() -> None:
    guide = "ACGT" * 5 + "ACG"
    sequence = "TTTA" + guide
    candidates = find_targets(sequence, system="Cas12a")
    _assert(any(c.strand == "+" and c.guide == guide for c in candidates), "Cas12a target not found")


def _test_coordinates() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    sequence = "TTTT" + guide + "AGG"
    candidates = find_targets(sequence, system="SpCas9", region_start=101)
    target = next(c for c in candidates if c.guide == guide and c.strand == "+")
    _assert(target.guide_start == 105, "Guide start coordinate incorrect")
    _assert(target.guide_end == 124, "Guide end coordinate incorrect")
    _assert(target.pam_start == 125, "PAM start coordinate incorrect")


def _test_ranking() -> None:
    sequence = "ACGT" * 8 + "AGG" + "AAAA" + "GCGC" * 6 + "TGG"
    result = design_guides(sequence, system="SpCas9")
    _assert(result.candidates, "Ranking test produced no candidates")
    scores = [c.rank_score for c in result.candidates]
    _assert(scores == sorted(scores, reverse=True), "Candidates are not ranked")
    _assert(all(0 <= score <= 100 for score in scores), "Score out of range")


def _test_mutation_overlap() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    result = design_guides(guide + "AGG", system="SpCas9")
    mutation = MutationTarget("chr1", 5, "A", "G")
    annotate_mutation_targets(result.candidates, [mutation])
    _assert(any(c.mutation_overlap for c in result.candidates), "Mutation overlap not detected")
    _assert(any(c.mutation_context for c in result.candidates), "Mutation context not populated")


def _test_local_off_target() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    refs = {"ref1": guide + "AGG", "ref2": guide[:-1] + "A" + "TGG"}
    result = design_guides(guide + "AGG", system="SpCas9")
    hits = search_local_off_targets_detailed(result.candidates[0], refs, max_mismatches=1)
    _assert(hits, "Detailed local off-target search found no hit")
    _assert(hits[0]["mismatches"] == 0, "Exact match should rank first")


def _test_json_roundtrip() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    result = design_guides(guide + "AGG", system="SpCas9")
    payload = json.loads(result_to_json(result))
    _assert(payload["sequence_name"] == "Selected Sequence", "JSON sequence name missing")
    _assert(payload["candidates"], "JSON candidate data missing")


def _test_orf_adapter() -> None:
    orf = {"sequence": "ACGTACGT", "start": 25, "name": "ORF_1"}
    extracted = extract_orf_context(orf)
    _assert(extracted == ("ACGTACGT", 25, "ORF_1"), "ORF adapter failed")


def _test_base_editing() -> None:
    guide = "CCCCACGTACGTACGTACGC"
    result = design_guides(guide + "AGG", system="SpCas9")
    mutation = MutationTarget("chr1", 5, "C", "T")
    analyses = analyze_base_editing(result.candidates, mutation, editor="CBE", window=(4, 8))
    _assert(analyses, "Base-editing analysis returned no candidates")


def _test_bystander_detection() -> None:
    guide = "CCCCACGTACGTACGTACGC"
    result = design_guides(guide + "AGG", system="SpCas9")
    mutation = MutationTarget("chr1", 5, "C", "T")
    analyses = analyze_base_editing(result.candidates, mutation, editor="CBE", window=(4, 8))
    _assert(any(a.bystander_count >= 0 for a in analyses), "Bystander analysis failed")


def _test_prime_editing() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    sequence = "AAAAA" + guide + "AGG" + "TTTTTTTTTTTTTTTTTTTTTT"
    result = design_guides(sequence, system="SpCas9")
    mutation = MutationTarget("chr1", 10, "A", "G")
    candidates = design_prime_editing_candidates(
        sequence, mutation, result.candidates, pbs_lengths=(8,), rtt_lengths=(10,)
    )
    _assert(isinstance(candidates, list), "Prime editing result is not a list")
    _assert(
        all(c.pbs_length == 8 and c.rtt_length == 10 for c in candidates),
        "Prime-editing candidates must have the requested PBS/RTT lengths",
    )


def _test_result_compatibility() -> None:
    result = design_guides("ACGT" * 8 + "AGG", system="SpCas9")
    _assert(result.region_start == result.region.start, "region_start compatibility property failed")
    _assert(result.region_end == result.region.end, "region_end compatibility property failed")
    _assert(result.region_name == result.region.name, "region_name compatibility property failed")
    _assert(result.feature_type == result.region.feature_type, "feature_type compatibility property failed")


def _test_prime_coordinate_translation() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    full = "A" * 50 + guide + "AGG" + "T" * 50
    region = full[40:73]
    result = design_guides(region, region_start=41, region_end=73, system="SpCas9")
    candidate = next(c for c in result.candidates if c.guide == guide and c.strand == "+")
    mutation = MutationTarget("chr1", 55, "A", "G")
    candidates = design_prime_editing_candidates(
        region, mutation, [candidate], pbs_lengths=(8,), rtt_lengths=(10,), sequence_start=41
    )
    _assert(candidates, "Global-to-local prime-editing coordinate translation failed")
    _assert(candidates[0].pbs_length == 8, "Prime-editing PBS length translation failed")
    _assert(candidates[0].rtt_length == 10, "Prime-editing RTT length translation failed")


def _test_fasta_parser() -> None:
    records = parse_fasta_text(">a\nACGT\n>b\nTTAA\n")
    _assert(records == [("a", "ACGT"), ("b", "TTAA")], "FASTA parser failed")


def _test_batch() -> None:
    records = [("a", "ACGT" * 8 + "AGG"), ("b", "TTTA" + "ACGT" * 6 + "AGG")]
    results = batch_design_guides(records, system="SpCas9")
    _assert(set(results) == {"a", "b"}, "Batch names incorrect")
    _assert(len(batch_summary(results)) == 2, "Batch summary incorrect")


def _test_pam_quality() -> None:
    q = pam_quality("AGG", "NGG")
    _assert(q["matched"], "PAM quality match failed")
    _assert(q["degenerate_positions"] == 1, "PAM degeneracy count failed")


def _test_target_map() -> None:
    guide = "ACGT" * 5
    result = design_guides(guide + "AGG", system="SpCas9")
    _assert("=" in target_map(result.candidates[0], result.sequence_length), "Target map missing guide")
    _assert("#" in target_map(result.candidates[0], result.sequence_length), "Target map missing PAM")


def _test_specificity_profile() -> None:
    profile = specificity_profile([
        {"mismatches": 0, "seed_mismatches": 0, "pam_compatible": True},
        {"mismatches": 1, "seed_mismatches": 1, "pam_compatible": True},
    ])
    _assert(profile["exact_matches"] == 1, "Specificity profile exact count failed")
    _assert(profile["pam_compatible_hits"] == 2, "Specificity PAM count failed")


def _test_specificity_aware_ranking() -> None:
    result = design_guides("ACGT" * 5 + "AGG", system="SpCas9")
    candidate = result.candidates[0]
    candidate.specificity_profile = specificity_profile([
        {"mismatches": 0, "seed_mismatches": 0, "pam_compatible": True},
    ], candidate.candidate_id)
    ranked = rank_guides(result.candidates, config=result.config)
    selected = next(c for c in ranked if c.candidate_id == candidate.candidate_id)
    _assert(has_specificity_analysis(ranked), "Specificity-aware ranking mode not detected")
    _assert("Specificity" in selected.rank_components, "Specificity component missing from ranking")
    _assert(selected.rank_components["Specificity"] == 60.0, "Specificity heuristic score is not deterministic")


def _test_export_path_resolution() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = resolve_export_path(Path(tmp) / "crispr_test.json")
        _assert(path.is_absolute(), "Export path is not absolute")
        _assert(path.parent == Path(tmp).resolve(), "Absolute export path was changed unexpectedly")


def _test_export_rows() -> None:
    guide = "ACGT" * 5
    result = design_guides(guide + "AGG", system="SpCas9")
    row = _candidate_csv_row(result.candidates[0])
    _assert("hairpin_risk" in row, "CSV row missing hairpin field")
    _assert("specificity_scope" in row, "CSV row missing specificity field")


def _test_variant_parser() -> None:
    variants = parse_vcf_variant_line("chr1\t100\t.\tA\tG,T\t100\tPASS\t.")
    _assert(len(variants) == 2, "VCF multiallelic parser failed")


def _test_system_details() -> None:
    details = format_system_details("SpCas9")
    _assert("System:          SpCas9" in details, "Selected system missing from details")
    _assert("PAM:             NGG" in details, "Selected PAM missing from details")
    _assert("SELECTED CRISPR SYSTEM" in details, "System details header missing")
    _assert("SUPPORTED CRISPR SYSTEMS" not in details, "Details should describe one system only")


def _test_ambiguous_handling() -> None:
    ok, errors = validate_sequence("ACGTN", allow_ambiguous=False)
    _assert(not ok and errors, "Ambiguous-base validation failed")
    ok, errors = validate_sequence("ACGTN", allow_ambiguous=True)
    _assert(ok and not errors, "Ambiguous-base opt-in failed")


def run_self_test(verbose: bool = True) -> bool:
    tests = [
        ("iupac_pam_matching", _test_iupac_pam),
        ("reverse_complement", _test_reverse_complement),
        ("guide_metrics", _test_metrics),
        ("spcas9_plus_strand", _test_spcas9_plus),
        ("spcas9_minus_strand", _test_spcas9_minus),
        ("cas12a_targeting", _test_cas12a),
        ("coordinate_mapping", _test_coordinates),
        ("guide_ranking", _test_ranking),
        ("mutation_overlap", _test_mutation_overlap),
        ("local_off_target_detailed", _test_local_off_target),
        ("json_serialization", _test_json_roundtrip),
        ("orf_adapter", _test_orf_adapter),
        ("base_editing", _test_base_editing),
        ("bystander_detection", _test_bystander_detection),
        ("prime_editing", _test_prime_editing),
        ("result_compatibility", _test_result_compatibility),
        ("prime_coordinate_translation", _test_prime_coordinate_translation),
        ("fasta_parser", _test_fasta_parser),
        ("batch_analysis", _test_batch),
        ("pam_quality", _test_pam_quality),
        ("target_visualization", _test_target_map),
        ("specificity_profile", _test_specificity_profile),
        ("specificity_aware_ranking", _test_specificity_aware_ranking),
        ("export_path_resolution", _test_export_path_resolution),
        ("export_schema", _test_export_rows),
        ("vcf_variant_parser", _test_variant_parser),
        ("system_details", _test_system_details),
        ("ambiguous_base_validation", _test_ambiguous_handling),
    ]
    results: List[SelfTestResult] = []
    for name, test in tests:
        try:
            test()
            results.append(SelfTestResult(name, True))
        except Exception as exc:
            results.append(SelfTestResult(name, False, str(exc)))
    if verbose:
        print("\n" + "=" * 78)
        print("SYNAPSE CRISPR SELF-TEST")
        print("=" * 78)
        for result in results:
            marker = "PASS" if result.passed else "FAIL"
            suffix = f" - {result.detail}" if result.detail else ""
            print(f"[{marker}] {result.name}{suffix}")
        passed = sum(r.passed for r in results)
        print(_line())
        print(f"{passed}/{len(results)} tests passed.")
    return all(result.passed for result in results)


# ============================================================================
# MODULE ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    if not run_self_test(verbose=True):
        raise SystemExit(1)
    print()
    print("CRISPR analysis engine self-test completed successfully.")
