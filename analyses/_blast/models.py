"""Data model of BLAST queries, parameters, runs, hits and databases.

``rank_hits`` lives here because ``BlastRun.ranked_hits`` uses it."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Sequence, Tuple

from analyses._blast.config import (
    DATABASES,
    DEFAULT_EMAIL,
    DEFAULT_MAX_HITS,
    DEFAULT_TOOL,
    PROGRAMS,
)


# ============================================================================
# DATA MODELS
# ============================================================================

@dataclass(frozen=True)
class BlastParameters:
    """Reproducible search configuration."""

    program: str = "blastn"
    database: str = "nt"
    evalue: float = 10.0
    max_target_seqs: int = DEFAULT_MAX_HITS
    word_size: Optional[int] = None
    gap_open: Optional[int] = None
    gap_extend: Optional[int] = None
    matrix: Optional[str] = None
    filter_query: Optional[str] = None
    num_threads: Optional[int] = None
    task: Optional[str] = None
    tool: str = DEFAULT_TOOL
    email: str = DEFAULT_EMAIL

    def normalized(self) -> "BlastParameters":
        program = self.program.lower().strip()
        database = self.database.strip()
        if program not in PROGRAMS:
            raise BlastAnalysisError(f"Unsupported BLAST program: {program}")
        if database not in DATABASES:
            raise BlastAnalysisError(f"Unsupported database: {database}")
        if self.evalue <= 0:
            raise BlastAnalysisError("E-value must be greater than zero.")
        if self.max_target_seqs < 1:
            raise BlastAnalysisError("max_target_seqs must be at least 1.")
        if self.num_threads is not None and self.num_threads < 1:
            raise BlastAnalysisError("num_threads must be at least 1.")

        expected_db_type = PROGRAMS[program]["database_type"]
        actual_db_type = DATABASES[database]["type"]
        if expected_db_type != actual_db_type:
            raise BlastAnalysisError(
                f"{program} requires a {expected_db_type} database, "
                f"but {database} is a {actual_db_type} database."
            )

        return BlastParameters(
            program=program,
            database=database,
            evalue=float(self.evalue),
            max_target_seqs=int(self.max_target_seqs),
            word_size=self.word_size,
            gap_open=self.gap_open,
            gap_extend=self.gap_extend,
            matrix=self.matrix,
            filter_query=self.filter_query,
            num_threads=self.num_threads,
            task=self.task.strip() if self.task else None,
            tool=self.tool.strip() or DEFAULT_TOOL,
            email=self.email.strip(),
        )


@dataclass
class BlastHSP:
    """One high-scoring pair (HSP) from a BLAST hit."""

    number: int
    bit_score: float
    raw_score: Optional[float]
    evalue: float
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    query_frame: Optional[int]
    subject_frame: Optional[int]
    identity: int
    positive: int
    gaps: int
    alignment_length: int
    query_sequence: str
    subject_sequence: str
    midline: str
    query_length: Optional[int] = None

    @property
    def identity_percent(self) -> float:
        if self.alignment_length <= 0:
            return 0.0
        return 100.0 * self.identity / self.alignment_length

    @property
    def positive_percent(self) -> float:
        if self.alignment_length <= 0:
            return 0.0
        return 100.0 * self.positive / self.alignment_length

    @property
    def gap_percent(self) -> float:
        if self.alignment_length <= 0:
            return 0.0
        return 100.0 * self.gaps / self.alignment_length

    @property
    def query_interval(self) -> Tuple[int, int]:
        """1-based inclusive query interval covered by this HSP."""
        low, high = sorted((self.query_start, self.query_end))
        return low, high

    @property
    def query_coverage(self) -> float:
        """Percent of the QUERY covered by this HSP (NCBI definition).

        Before 2026-10-29 this property returned the non-gap fraction of the
        alignment, which is 100 % for every gap-free HSP regardless of how
        much of the query aligns (audit P0-1). That value is still available
        as ``alignment_nongap_percent``. Returns NaN when the query length is
        unknown.
        """
        if not self.query_length:
            return math.nan
        low, high = self.query_interval
        if low < 1:
            return math.nan
        covered = high - low + 1
        return min(100.0, 100.0 * covered / self.query_length)

    @property
    def alignment_nongap_percent(self) -> float:
        """Percent of aligned query columns that are not gaps."""
        if not self.query_sequence:
            return 0.0
        non_gap = sum(base != "-" for base in self.query_sequence)
        return 100.0 * non_gap / max(1, len(self.query_sequence))


@dataclass
class BlastHit:
    """A BLAST subject hit containing one or more HSPs."""

    hit_id: str
    accession: str
    definition: str
    length: Optional[int]
    hsps: List[BlastHSP] = field(default_factory=list)
    scientific_names: List[str] = field(default_factory=list)
    taxonomy_ids: List[str] = field(default_factory=list)
    query_length: Optional[int] = None

    def __post_init__(self) -> None:
        if self.query_length:
            for hsp in self.hsps:
                if not hsp.query_length:
                    hsp.query_length = self.query_length

    @property
    def combined_query_coverage_percent(self) -> float:
        """Percent of the query covered by the union of all HSPs (NCBI "Query cover")."""
        length = self.query_length or next(
            (h.query_length for h in self.hsps if h.query_length), None
        )
        if not length or not self.hsps:
            return math.nan if not length else 0.0
        intervals = sorted(h.query_interval for h in self.hsps)
        covered = 0
        current_low, current_high = intervals[0]
        for low, high in intervals[1:]:
            if low <= current_high + 1:
                current_high = max(current_high, high)
            else:
                covered += current_high - current_low + 1
                current_low, current_high = low, high
        covered += current_high - current_low + 1
        return min(100.0, 100.0 * covered / length)

    @property
    def best_hsp(self) -> Optional[BlastHSP]:
        if not self.hsps:
            return None
        return max(self.hsps, key=lambda h: (h.bit_score, -h.evalue))

    @property
    def best_bit_score(self) -> float:
        return self.best_hsp.bit_score if self.best_hsp else 0.0

    @property
    def best_evalue(self) -> Optional[float]:
        return self.best_hsp.evalue if self.best_hsp else None

    @property
    def best_identity_percent(self) -> float:
        return self.best_hsp.identity_percent if self.best_hsp else 0.0

    @property
    def best_coverage_percent(self) -> float:
        """Query coverage of the hit: union of all HSPs over the query length.

        This matches NCBI's "Query cover". The best single HSP's coverage is
        ``best_hsp_coverage_percent``.
        """
        return self.combined_query_coverage_percent if self.hsps else 0.0

    @property
    def best_hsp_coverage_percent(self) -> float:
        return self.best_hsp.query_coverage if self.best_hsp else 0.0

    @property
    def total_alignment_length(self) -> int:
        return sum(h.alignment_length for h in self.hsps)


@dataclass
class BlastQuery:
    """Query sequence and provenance information."""

    name: str
    sequence: str
    sequence_type: str = "DNA"
    source: str = "Selected Sequence"
    region_start: int = 1
    region_end: Optional[int] = None
    feature_type: str = "sequence"
    feature_name: Optional[str] = None
    notes: List[str] = field(default_factory=list)
    # Set only when the sequence itself is unavailable, as for a run resumed
    # from its request ID: the BLAST XML reports the query length but not the
    # query sequence (P1-7).
    reported_length: Optional[int] = None

    @property
    def length(self) -> int:
        if self.sequence:
            return len(self.sequence)
        return int(self.reported_length or 0)


@dataclass
class BlastRunMetadata:
    """Metadata required to make a remote BLAST run reproducible."""

    rid: str
    submitted_at_utc: str
    completed_at_utc: Optional[str]
    program: str
    database: str
    database_type: str
    query_name: str
    query_length: int
    query_source: str
    parameters: Dict[str, Any]
    ncbi_status: str
    ncbi_version: Optional[str] = None
    ncbi_database_reported: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    engine: str = "NCBI Remote"
    local_database_path: Optional[str] = None
    # Effective parameters as reported by BLAST in its XML output
    # (expect, scoring, gap costs, filter) - the values actually used.
    reported_parameters: Dict[str, Any] = field(default_factory=dict)
    # Database provenance (local: blastdbcmd -info; remote: NCBI statistics).
    database_details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class BlastRun:
    """Complete Synapse BLAST result."""

    metadata: BlastRunMetadata
    query: BlastQuery
    hits: List[BlastHit] = field(default_factory=list)
    database_statistics: Dict[str, Any] = field(default_factory=dict)
    interpretation_notes: List[str] = field(default_factory=list)

    @property
    def hit_count(self) -> int:
        return len(self.hits)

    @property
    def significant_hit_count(self) -> int:
        return sum(
            1 for hit in self.hits
            if hit.best_evalue is not None
            and hit.best_evalue <= self.metadata.parameters.get("evalue", 10.0)
        )

    @property
    def best_hit(self) -> Optional[BlastHit]:
        return rank_hits(self.hits)[0] if self.hits else None


@dataclass
class HitEvidence:
    """Transparent post-BLAST evidence profile."""

    hit_id: str
    accession: str
    identity_percent: float
    query_coverage_percent: float
    evalue: Optional[float]
    bit_score: float
    alignment_length: int
    evidence_score: float
    evidence_band: str
    notes: List[str] = field(default_factory=list)


@dataclass
class BlastComparison:
    """Comparison of two independent BLAST runs."""

    query_a: str
    query_b: str
    common_accessions: List[str]
    unique_to_a: List[str]
    unique_to_b: List[str]
    identity_differences: Dict[str, float]
    coverage_differences: Dict[str, float]


@dataclass
class SelfTestResult:
    name: str
    passed: bool
    detail: str = ""


class BlastAnalysisError(Exception):
    """Raised for invalid input, API errors, or malformed BLAST results."""


# ============================================================================
# GENERAL HELPERS
# ============================================================================

def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _evalue_text(evalue: Optional[float]) -> str:
    """An E-value with three significant figures, or n/a when absent."""
    return "n/a" if evalue is None else f"{evalue:.3g}"


# ============================================================================
# LOCAL BLAST+ INFRASTRUCTURE
# ============================================================================

@dataclass
class LocalBlastDatabase:
    """Detected local BLAST database and verification metadata."""

    name: str
    database_type: str
    path: str
    source: str = "Unknown"
    installed: bool = False
    structurally_valid: bool = False
    verification_status: str = "NOT_CHECKED"
    verification_message: str = ""
    file_count: int = 0
    size_bytes: int = 0
    modified_at: Optional[str] = None
    title: Optional[str] = None
    sequence_count: Optional[int] = None
    total_length: Optional[int] = None
    blastdb_version: Optional[str] = None
    database_date: Optional[str] = None
    volumes: List[str] = field(default_factory=list)
    volume_paths: List[str] = field(default_factory=list)
    volume_count: int = 0
    expected_volume_count: Optional[int] = None
    missing_volumes: List[str] = field(default_factory=list)
    incomplete_volumes: List[str] = field(default_factory=list)
    missing_volume_files: Dict[str, List[str]] = field(default_factory=dict)
    manifest_metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def size_mb(self) -> float:
        return self.size_bytes / (1024 * 1024)

    @property
    def size_gb(self) -> float:
        return self.size_bytes / (1024 * 1024 * 1024)

    @property
    def exists(self) -> bool:
        return self.installed and self.structurally_valid


def rank_hits(hits: Sequence[BlastHit]) -> List[BlastHit]:
    """
    Deterministic ranking for display.

    Primary criteria:
    1. bit score descending
    2. E-value ascending
    3. identity descending
    4. query coverage descending
    """
    return sorted(
        hits,
        key=lambda hit: (
            -hit.best_bit_score,
            hit.best_evalue if hit.best_evalue is not None else math.inf,
            -hit.best_identity_percent,
            -hit.best_coverage_percent,
            hit.accession or hit.hit_id,
        ),
    )
