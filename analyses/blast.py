"""
Synapse - BLAST Search Engine
=============================

Professional BLAST integration and result-analysis module for Project Synapse.

Scientific design
-----------------
This module does NOT reimplement the BLAST algorithm and does not fabricate
database hits. It uses the NCBI BLAST URL API for remote searches and the official
NCBI BLAST+ command-line applications for local searches. Synapse performs
transparent post-processing and database management around those engines.

Supported BLAST programs
------------------------
- blastn   : nucleotide query -> nucleotide database
- blastp   : protein query   -> protein database
- blastx   : translated nucleotide query -> protein database
- tblastn  : protein query -> translated nucleotide database
- tblastx  : translated nucleotide query -> translated nucleotide database

Core capabilities
-----------------
- NCBI remote BLAST submission and polling.
- Local NCBI BLAST+ execution against managed databases.
- Local BLAST database status, verification, download/update and removal.
- Custom FASTA -> local BLAST database creation through makeblastdb.
- XML result parsing into typed Synapse records.
- Identity, coverage, E-value, bit score and alignment metrics.
- HSP-level alignment inspection.
- Transparent hit ranking / evidence profiling.
- Taxonomic grouping when taxonomy information is returned by BLAST.
- Reproducible run metadata and parameter recording.
- ORF-aware query preparation.
- Comparison of independent BLAST runs.
- JSON, CSV and text report export.
- Deterministic self-test suite that never calls the network.

Scientific boundaries
---------------------
- A BLAST hit is a sequence-similarity result, not proof of biological identity,
  function, pathogenicity, clinical relevance or evolutionary relationship.
- Interpretation is deliberately evidence-based and transparent.
- NCBI database contents and versions change over time; run metadata therefore
  records the requested database, program, parameters and timestamp.
- Remote BLAST requires internet access and is subject to NCBI service limits.
- Local BLAST requires an installed NCBI BLAST+ package and an installed local
  database.
- Official NCBI database downloads are explicit user actions; Synapse never
  downloads multi-gigabyte databases automatically.
- This module does not claim genome-wide biological conclusions beyond the
  searched database and returned alignments.

No third-party Python packages are required.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


# ============================================================================
# CONSTANTS
# ============================================================================

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
DEFAULT_TOOL = os.environ.get("SYNAPSE_NCBI_TOOL", "Synapse").strip() or "Synapse"
DEFAULT_EMAIL = os.environ.get("SYNAPSE_NCBI_EMAIL", "").strip()
API_TIMEOUT_SECONDS = 30
# NCBI asks API clients not to contact the BLAST server more often than once
# every 10 seconds and not to poll an individual RID more often than once a
# minute. These limits are enforced by the client instead of relying on callers.
NCBI_MIN_REQUEST_INTERVAL_SECONDS = 10.0
NCBI_MIN_POLL_INTERVAL_SECONDS = 60.0
DEFAULT_POLL_INTERVAL_SECONDS = NCBI_MIN_POLL_INTERVAL_SECONDS
DEFAULT_MAX_WAIT_SECONDS = 300.0
DEFAULT_MAX_HITS = 50


# Local BLAST database storage. The path can be moved without changing Synapse
# code by setting SYNAPSE_BLAST_DB_ROOT in the environment.
LOCAL_DB_ROOT_ENV = "SYNAPSE_BLAST_DB_ROOT"
LOCAL_DB_MANIFEST_NAME = "manifest.json"


def _resolve_default_local_db_root() -> Path:
    """Resolve Synapse's BLAST database root independently of the CWD.

    Synapse may import ``blast.py`` from a subdirectory such as ``analyses``.
    Using ``Path(__file__).parent / data / blast_db`` in that situation creates
    the wrong location (for example ``analyses/data/blast_db``). We therefore
    walk upward and prefer an existing project-level ``data/blast_db``.

    The environment variable remains the highest-priority override.
    """
    configured = os.environ.get(LOCAL_DB_ROOT_ENV, "").strip()
    if configured:
        return Path(configured).expanduser().resolve()

    module_dir = Path(__file__).resolve().parent
    candidates = [module_dir, *module_dir.parents]

    # Prefer an already existing BLAST root. This preserves databases that
    # were created before Synapse starts and avoids creating a second root.
    for candidate in candidates:
        blast_root = candidate / "data" / "blast_db"
        if blast_root.is_dir():
            return blast_root.resolve()

    # If no BLAST directory exists yet, identify the project root by common
    # Synapse markers before falling back to the module directory.
    for candidate in candidates:
        if (candidate / "main.py").is_file() or (candidate / "data").is_dir():
            return (candidate / "data" / "blast_db").resolve()

    return (module_dir / "data" / "blast_db").resolve()


DEFAULT_LOCAL_DB_ROOT = _resolve_default_local_db_root()

# NCBI preformatted database catalog. These names are passed directly to
# update_blastdb.pl. The catalog is deliberately conservative: Synapse never
# downloads a large database automatically.
LOCAL_DATABASE_CATALOG = {
    "refseq_protein": {
        "type": "protein",
        "description": "NCBI RefSeq Protein",
        "source": "NCBI preformatted BLAST database",
        "recommended": True,
    },
    "refseq_rna": {
        "type": "nucleotide",
        "description": "NCBI RefSeq RNA",
        "source": "NCBI preformatted BLAST database",
        "recommended": True,
    },
    "refseq_genomic": {
        "type": "nucleotide",
        "description": "NCBI RefSeq genomic",
        "source": "NCBI preformatted BLAST database",
        "recommended": False,
    },
    "swissprot": {
        "type": "protein",
        "description": "UniProtKB/Swiss-Prot via NCBI",
        "source": "NCBI preformatted BLAST database",
        "recommended": True,
    },
    "nt": {
        "type": "nucleotide",
        "description": "NCBI nucleotide collection",
        "source": "NCBI preformatted BLAST database",
        "recommended": False,
    },
    "nr": {
        "type": "protein",
        "description": "NCBI non-redundant protein collection",
        "source": "NCBI preformatted BLAST database",
        "recommended": False,
    },
}

BLAST_EXECUTABLES = (
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "blastdbcmd",
    "makeblastdb",
)
BLAST_UPDATE_SCRIPT = "update_blastdb.pl"

DNA_BASES = set("ACGT")
PROTEIN_AMINO_ACIDS = set("ABCDEFGHIKLMNPQRSTVWXYZ*")
AMBIGUOUS_DNA = set("NRYKMSWBDHV")

PROGRAMS = {
    "blastn": {
        "query_type": "DNA",
        "database_type": "nucleotide",
        "description": "Nucleotide query against a nucleotide database.",
    },
    "blastp": {
        "query_type": "protein",
        "database_type": "protein",
        "description": "Protein query against a protein database.",
    },
    "blastx": {
        "query_type": "DNA",
        "database_type": "protein",
        "description": "Six-frame translated nucleotide query against protein.",
    },
    "tblastn": {
        "query_type": "protein",
        "database_type": "nucleotide",
        "description": "Protein query against six-frame translated nucleotide sequences.",
    },
    "tblastx": {
        "query_type": "DNA",
        "database_type": "nucleotide",
        "description": "Translated nucleotide query against translated nucleotide sequences.",
    },
}

DATABASES = {
    "nt": {
        "type": "nucleotide",
        "description": "NCBI nucleotide collection.",
    },
    "refseq_rna": {
        "type": "nucleotide",
        "description": "NCBI RefSeq RNA collection.",
    },
    "refseq_genomic": {
        "type": "nucleotide",
        "description": "NCBI RefSeq genomic collection.",
    },
    "nr": {
        "type": "protein",
        "description": "NCBI non-redundant protein collection.",
    },
    "refseq_protein": {
        "type": "protein",
        "description": "NCBI RefSeq protein collection.",
    },
    "swissprot": {
        "type": "protein",
        "description": "UniProtKB/Swiss-Prot protein collection as exposed by NCBI.",
    },
}

# NCBI BLAST XML 2/legacy XML namespaces can vary. We deliberately parse
# local tag names so the parser is resilient to namespace prefixes.
BLAST_TAG_ALIASES = {
    "Iteration",
    "Hit",
    "Hsp",
    "Hit_id",
    "Hit_def",
    "Hit_accession",
    "Hit_len",
    "Hsp_num",
    "Hsp_bit-score",
    "Hsp_score",
    "Hsp_evalue",
    "Hsp_query-from",
    "Hsp_query-to",
    "Hsp_hit-from",
    "Hsp_hit-to",
    "Hsp_query-frame",
    "Hsp_hit-frame",
    "Hsp_identity",
    "Hsp_positive",
    "Hsp_gaps",
    "Hsp_align-len",
    "Hsp_qseq",
    "Hsp_hseq",
    "Hsp_midline",
    "BlastOutput_program",
    "BlastOutput_version",
    "BlastOutput_db",
    "Parameters_expect",
    "Parameters_matrix",
    "Parameters_gap-open",
    "Parameters_gap-extend",
    "Statistics_db-num",
    "Statistics_db-len",
    "Statistics_hsp-len",
}


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
    def query_coverage(self) -> float:
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

    @property
    def length(self) -> int:
        return len(self.sequence)


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


def _local_tag(element: ET.Element) -> str:
    """Return XML tag without namespace."""
    return element.tag.rsplit("}", 1)[-1]


def _child_text(element: ET.Element, tag: str, default: str = "") -> str:
    for child in element:
        if _local_tag(child) == tag:
            return (child.text or default).strip()
    return default


def _descendants(element: ET.Element, tag: str) -> Iterable[ET.Element]:
    for node in element.iter():
        if _local_tag(node) == tag:
            yield node


def _to_int(value: str, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _to_float(value: str, default: Optional[float] = None) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def normalize_sequence(sequence: str, sequence_type: str = "DNA") -> str:
    """Normalize a sequence while preserving biological alphabet semantics."""
    if not isinstance(sequence, str):
        raise BlastAnalysisError("Sequence must be a string.")

    cleaned = "".join(sequence.split()).upper()
    if not cleaned:
        raise BlastAnalysisError("Sequence is empty.")

    sequence_type = sequence_type.upper()
    if sequence_type == "DNA":
        allowed = DNA_BASES | AMBIGUOUS_DNA
        invalid = sorted(set(cleaned) - allowed)
        if invalid:
            raise BlastAnalysisError(
                f"Invalid DNA characters: {', '.join(invalid)}"
            )
    elif sequence_type == "PROTEIN":
        invalid = sorted(set(cleaned) - PROTEIN_AMINO_ACIDS)
        if invalid:
            raise BlastAnalysisError(
                f"Invalid protein characters: {', '.join(invalid)}"
            )
    else:
        raise BlastAnalysisError("sequence_type must be DNA or PROTEIN.")

    return cleaned


def infer_sequence_type(sequence: str) -> str:
    """Infer DNA vs protein conservatively."""
    cleaned = "".join(sequence.split()).upper()
    if not cleaned:
        raise BlastAnalysisError("Cannot infer type from an empty sequence.")
    if set(cleaned) <= (DNA_BASES | AMBIGUOUS_DNA):
        return "DNA"
    return "PROTEIN"


def build_query(
    sequence: str,
    *,
    name: str = "Synapse_Query",
    sequence_type: Optional[str] = None,
    source: str = "Selected Sequence",
    region_start: int = 1,
    region_end: Optional[int] = None,
    feature_type: str = "sequence",
    feature_name: Optional[str] = None,
) -> BlastQuery:
    inferred = sequence_type.upper() if sequence_type else infer_sequence_type(sequence)
    normalized = normalize_sequence(sequence, inferred)
    end = region_end if region_end is not None else region_start + len(normalized) - 1
    if region_start < 1 or end < region_start:
        raise BlastAnalysisError("Invalid query coordinates.")

    return BlastQuery(
        name=name.strip() or "Synapse_Query",
        sequence=normalized,
        sequence_type=inferred,
        source=source,
        region_start=region_start,
        region_end=end,
        feature_type=feature_type,
        feature_name=feature_name,
    )


def _resolve_synapse_user_path(path: str | Path) -> Path:
    """Resolve a user-supplied path independently of the process CWD.

    Synapse can be launched from Windows Terminal, VS Code, a shortcut, or
    another application. In some launch modes Windows uses ``C:\\Windows\\System32``
    as the current working directory. Relative paths entered by the user
    should nevertheless be interpreted relative to the Synapse project root.

    Absolute paths are preserved. Relative paths are first interpreted from
    the Synapse project root (the directory containing ``analyses``).
    """
    candidate = Path(path).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()

    working_directory = Path.cwd()
    # Windows shortcuts can launch console applications with System32 as the
    # process CWD. In that special case, use Synapse's application root as a
    # practical fallback; otherwise preserve normal relative-path semantics.
    if os.name == "nt" and working_directory.resolve() == Path(os.environ.get("SystemRoot", r"C:\Windows")).resolve() / "System32":
        working_directory = Path(__file__).resolve().parent.parent
    return (working_directory / candidate).resolve()


def load_fasta_queries(
    fasta_path: str | Path,
    *,
    source: str = "FASTA File",
) -> List[BlastQuery]:
    """Load one or more BLAST queries from a FASTA file.

    The parser is intentionally small and dependency-free because Synapse only
    needs FASTA headers and sequence records for BLAST query preparation.
    Every record is validated through :func:`build_query`, so invalid DNA or
    protein characters are rejected before a BLAST search is started.

    A multi-record FASTA is returned as multiple ``BlastQuery`` objects. The
    interactive menu lets the user explicitly choose which record to search,
    avoiding an implicit first-record-only behavior.
    """
    path = _resolve_synapse_user_path(fasta_path)
    if not path.is_file():
        raise BlastAnalysisError(f"FASTA file not found: {path}")

    records: List[Tuple[str, str]] = []
    current_header: Optional[str] = None
    current_sequence: List[str] = []

    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except OSError as exc:
        raise BlastAnalysisError(f"Could not read FASTA file: {exc}") from exc
    except UnicodeError as exc:
        raise BlastAnalysisError(f"FASTA file is not valid UTF-8: {exc}") from exc

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(current_sequence)))
            current_header = line[1:].strip() or f"FASTA_Record_{len(records) + 1}"
            current_sequence = []
        else:
            if current_header is None:
                raise BlastAnalysisError(
                    "Invalid FASTA file: sequence data appears before the first header."
                )
            current_sequence.append(line)

    if current_header is not None:
        records.append((current_header, "".join(current_sequence)))

    if not records:
        raise BlastAnalysisError("FASTA file contains no sequence records.")

    queries: List[BlastQuery] = []
    for index, (header, sequence) in enumerate(records, start=1):
        if not sequence:
            raise BlastAnalysisError(
                f"FASTA record {index} ('{header}') contains no sequence."
            )
        queries.append(
            build_query(
                sequence,
                name=header,
                source=source,
            )
        )

    return queries


# ============================================================================
# NCBI BLAST API
# ============================================================================

class NCBIBlastClient:
    """
    Minimal dependency-free client for the NCBI BLAST URL API.

    The client performs:
        PUT -> receive RID
        GET -> poll status
        GET -> retrieve XML

    It intentionally does not silently fall back to fake/local hits.
    """

    def __init__(
        self,
        *,
        tool: str = DEFAULT_TOOL,
        email: str = DEFAULT_EMAIL,
        timeout: float = API_TIMEOUT_SECONDS,
        user_agent: Optional[str] = None,
        min_request_interval: float = NCBI_MIN_REQUEST_INTERVAL_SECONDS,
        min_poll_interval: float = NCBI_MIN_POLL_INTERVAL_SECONDS,
    ) -> None:
        raw_tool = tool.strip() or DEFAULT_TOOL
        # NCBI recommends a tool identifier without internal spaces.
        self.tool = re.sub(r"\s+", "_", raw_tool)
        self.email = email.strip()
        self.timeout = float(timeout)
        self.min_request_interval = max(0.0, float(min_request_interval))
        self.min_poll_interval = max(1.0, float(min_poll_interval))
        self.user_agent = user_agent or f"{self.tool}/1.0"
        self._last_request_at: Optional[float] = None

    def _throttle_request(self) -> None:
        """Respect NCBI's public BLAST request-frequency guidance."""
        if self._last_request_at is None:
            return
        elapsed = time.monotonic() - self._last_request_at
        remaining = self.min_request_interval - elapsed
        if remaining > 0:
            time.sleep(remaining)

    def _request(
        self,
        params: Mapping[str, Any],
        *,
        method: str = "GET",
    ) -> str:
        encoded = urlencode(
            [(key, value) for key, value in params.items() if value is not None]
        )

        if method.upper() == "POST":
            request = Request(
                NCBI_BLAST_URL,
                data=encoded.encode("utf-8"),
                method="POST",
            )
        else:
            separator = "&" if "?" in NCBI_BLAST_URL else "?"
            request = Request(
                f"{NCBI_BLAST_URL}{separator}{encoded}",
                method="GET",
            )

        request.add_header("User-Agent", self.user_agent)
        request.add_header("Accept", "text/plain, application/xml, */*")

        self._throttle_request()
        try:
            with urlopen(request, timeout=self.timeout) as response:
                self._last_request_at = time.monotonic()
                return response.read().decode("utf-8", errors="replace")
        except HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            raise BlastAnalysisError(
                f"NCBI BLAST HTTP error {exc.code}: {body[:500]}"
            ) from exc
        except URLError as exc:
            raise BlastAnalysisError(f"NCBI BLAST connection error: {exc}") from exc
        except TimeoutError as exc:
            raise BlastAnalysisError("NCBI BLAST request timed out.") from exc

    def submit(
        self,
        query: BlastQuery,
        parameters: BlastParameters,
    ) -> Tuple[str, Optional[float], str]:
        params = parameters.normalized()
        payload: Dict[str, Any] = {
            "CMD": "Put",
            "PROGRAM": params.program,
            "DATABASE": params.database,
            "QUERY": f">{query.name}\r\n{query.sequence}",
            "EXPECT": params.evalue,
            "HITLIST_SIZE": params.max_target_seqs,
            "TOOL": params.tool,
        }

        if params.email:
            payload["EMAIL"] = params.email
        if params.word_size is not None:
            payload["WORD_SIZE"] = params.word_size
        if params.gap_open is not None:
            payload["GAPOPEN"] = params.gap_open
        if params.gap_extend is not None:
            payload["GAPEXT"] = params.gap_extend
        if params.matrix:
            payload["MATRIX_NAME"] = params.matrix
        if params.filter_query:
            payload["FILTER"] = params.filter_query

        response = self._request(payload, method="POST")

        rid_match = re.search(r"RID\s*=\s*([A-Z0-9_-]+)", response, re.IGNORECASE)
        if not rid_match:
            raise BlastAnalysisError(
                "NCBI did not return a request ID (RID). "
                f"Response: {response[:1000]}"
            )

        rtoe_match = re.search(r"RTOE\s*=\s*(\d+)", response, re.IGNORECASE)
        rtoe = float(rtoe_match.group(1)) if rtoe_match else None
        rid = rid_match.group(1)
        return rid, rtoe, response

    def status(self, rid: str) -> Tuple[str, str]:
        response = self._request(
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo",
            }
        )

        status_match = re.search(
            r"Status\s*=\s*([A-Za-z]+)", response, re.IGNORECASE
        )
        status = status_match.group(1).upper() if status_match else "UNKNOWN"
        return status, response

    def fetch_xml(self, rid: str) -> str:
        response = self._request(
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "XML",
            }
        )
        if "<BlastOutput" not in response and "BlastOutput2" not in response:
            raise BlastAnalysisError(
                "NCBI returned no recognizable BLAST XML result. "
                f"Response: {response[:1000]}"
            )
        return response

    def run(
        self,
        query: BlastQuery,
        parameters: BlastParameters,
        *,
        poll_interval: float = DEFAULT_POLL_INTERVAL_SECONDS,
        max_wait: float = DEFAULT_MAX_WAIT_SECONDS,
    ) -> BlastRun:
        params = parameters.normalized()

        if params.email == "":
            # NCBI can accept requests without email, but Synapse exposes this
            # explicitly instead of inventing one.
            email_warning = (
                "No contact email was configured for the NCBI request. "
                "Set an email address in Synapse configuration when appropriate."
            )
        else:
            email_warning = None

        submitted_at = _utc_now()
        rid, rtoe, _submit_response = self.submit(query, params)

        # NCBI recommends no more than one status request per RID per minute.
        # RTOE is only an estimate, so we honor it while enforcing the minimum
        # polling interval as the hard lower bound.
        wait_time = max(0.0, rtoe or 0.0, self.min_poll_interval)
        if wait_time:
            time.sleep(min(wait_time, max_wait))

        start = time.monotonic()
        status = "WAITING"

        while True:
            status, _status_response = self.status(rid)

            if status == "READY":
                break
            if status in {"FAILED", "UNKNOWN"}:
                raise BlastAnalysisError(
                    f"NCBI BLAST request {rid} ended with status {status}."
                )
            if time.monotonic() - start >= max_wait:
                raise BlastAnalysisError(
                    f"BLAST request {rid} exceeded the maximum wait time "
                    f"of {max_wait:.0f} seconds."
                )

            time.sleep(max(self.min_poll_interval, poll_interval))

        xml_text = self.fetch_xml(rid)
        completed_at = _utc_now()

        run = parse_blast_xml(
            xml_text,
            query=query,
            rid=rid,
            parameters=params,
            submitted_at=submitted_at,
            completed_at=completed_at,
        )

        if email_warning:
            run.metadata.warnings.append(email_warning)

        return run



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


class BlastExecutableManager:
    """Locate BLAST+ executables without hard-coding a user-specific install path."""

    def __init__(self, blast_bin: Optional[str | Path] = None) -> None:
        self.explicit_bin = Path(blast_bin).expanduser() if blast_bin else None

    def _candidate_bins(self) -> List[Path]:
        candidates: List[Path] = []

        if self.explicit_bin:
            candidates.append(self.explicit_bin)

        for env_name in ("SYNAPSE_BLAST_BIN", "NCBI_BLAST_BIN"):
            value = os.environ.get(env_name, "").strip()
            if value:
                candidates.append(Path(value).expanduser())

        for env_name in ("BLAST_HOME", "NCBI_BLAST_HOME"):
            value = os.environ.get(env_name, "").strip()
            if value:
                root = Path(value).expanduser()
                candidates.extend([root / "bin", root])

        # Common Windows installation roots. These are only discovery hints;
        # no user-specific path is assumed.
        if os.name == "nt":
            for env_name in ("ProgramFiles", "ProgramFiles(x86)", "LOCALAPPDATA"):
                value = os.environ.get(env_name, "").strip()
                if not value:
                    continue
                base = Path(value)
                candidates.extend(
                    [
                        base / "NCBI" / "blast",
                        base / "NCBI",
                        base / "blast",
                    ]
                )
                try:
                    candidates.extend(
                        [path / "bin" for path in base.glob("NCBI/blast-*/bin")]
                    )
                except OSError:
                    pass

        # The executable may live next to this module in a portable Synapse
        # deployment or inside a local blast*/bin directory.
        module_root = Path(__file__).resolve().parent
        candidates.extend(
            [
                module_root / "bin",
                module_root / "blast" / "bin",
            ]
        )
        try:
            candidates.extend([path / "bin" for path in module_root.glob("blast-*/bin")])
        except OSError:
            pass

        unique: List[Path] = []
        seen: set[str] = set()
        for candidate in candidates:
            key = str(candidate.resolve()) if candidate.exists() else str(candidate)
            if key not in seen:
                seen.add(key)
                unique.append(candidate)
        return unique

    def resolve(self, executable: str) -> Optional[str]:
        """Return a usable executable path, or None if unavailable."""
        direct = shutil.which(executable)
        if direct:
            return direct

        names = [executable]
        if os.name == "nt" and not executable.lower().endswith(".exe"):
            names.insert(0, executable + ".exe")

        for directory in self._candidate_bins():
            for name in names:
                candidate = directory / name
                if candidate.is_file():
                    return str(candidate)

        return None

    def available(self) -> Dict[str, Optional[str]]:
        return {name: self.resolve(name) for name in BLAST_EXECUTABLES}

    def require(self, executable: str) -> str:
        path = self.resolve(executable)
        if not path:
            raise BlastAnalysisError(
                f"BLAST+ executable '{executable}' was not found. "
                "Install NCBI BLAST+ or configure SYNAPSE_BLAST_BIN."
            )
        return path

    def resolve_perl(self) -> Optional[str]:
        return shutil.which("perl")

    def resolve_update_script(self) -> Optional[str]:
        """Locate NCBI's update_blastdb.pl without assuming a user-specific path.

        On Windows, ``blastn.exe`` may be discoverable through PATH while the
        accompanying Perl script is not itself on PATH.  In that common NCBI
        installation layout, the most reliable discovery method is to resolve
        a BLAST+ executable first and inspect its sibling ``bin`` directory.
        """
        direct = shutil.which(BLAST_UPDATE_SCRIPT)
        if direct:
            return direct

        # Prefer the directory containing a BLAST+ executable already visible
        # to the operating system.  This handles standard Windows installs
        # such as ``...\\NCBI\\blast-2.17.0+\\bin`` even when the Perl
        # script itself is not included in PATH.
        for executable in BLAST_EXECUTABLES:
            resolved_executable = shutil.which(executable)
            if not resolved_executable:
                continue
            candidate = Path(resolved_executable).resolve().parent / BLAST_UPDATE_SCRIPT
            if candidate.is_file():
                return str(candidate)

        for directory in self._candidate_bins():
            candidate = directory / BLAST_UPDATE_SCRIPT
            if candidate.is_file():
                return str(candidate)

        return None


def _parse_blastdbcmd_info(output: str) -> Dict[str, Any]:
    """Extract stable summary fields from blastdbcmd -info output."""
    text = output or ""
    info: Dict[str, Any] = {
        "title": None,
        "sequence_count": None,
        "total_length": None,
        "blastdb_version": None,
        "volumes": [],
    }

    database_match = re.search(r"^Database:\s*(.+)$", text, re.MULTILINE)
    if database_match:
        info["title"] = database_match.group(1).strip()

    count_match = re.search(
        r"([\d,]+)\s+sequences;\s+([\d,]+)\s+total\s+bases",
        text,
        re.IGNORECASE,
    )
    if count_match:
        info["sequence_count"] = int(count_match.group(1).replace(",", ""))
        info["total_length"] = int(count_match.group(2).replace(",", ""))

    version_match = re.search(
        r"BLASTDB\s+Version:\s*([0-9]+)",
        text,
        re.IGNORECASE,
    )
    if version_match:
        info["blastdb_version"] = version_match.group(1)

    in_volumes = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.lower() == "volumes:":
            in_volumes = True
            continue
        if in_volumes and stripped:
            if stripped.startswith("Database:") or stripped.startswith("Date:"):
                break
            info["volumes"].append(stripped)

    return info


def _extract_blast_volume_index(volume_name: str) -> Optional[int]:
    """Return the numeric BLAST volume index while preserving zero padding elsewhere.

    NCBI alias files and extracted database layouts can represent the same
    volume as ``refseq_rna.0`` or ``refseq_rna.00``. The numeric index is the
    stable identity; the actual filesystem name (including zero padding) is
    preserved by the resolver.
    """
    match = re.search(r"\.(\d+)$", str(volume_name).strip())
    if not match:
        return None
    try:
        return int(match.group(1))
    except ValueError:
        return None


def _parse_blast_alias_volume_names(alias_path: Path) -> List[str]:
    """Read volume names declared by an NCBI BLAST alias (.nal/.pal) file.

    NCBI multi-volume databases can be stored in dedicated folders such as
    ``refseq_rna.00/`` while the alias file lives one directory above them.
    The alias remains the authoritative source for which volumes belong to the
    logical database, so Synapse uses DBLIST when it is available.
    """
    if not alias_path.is_file():
        return []

    try:
        payload = alias_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []

    match = re.search(
        r"(?ims)^\s*DBLIST\s+(.+?)(?=^\s*[A-Z][A-Z0-9_]*\s+|\Z)",
        payload,
    )
    if not match:
        return []

    block = match.group(1)
    values = re.findall(r'"([^"]+)"', block)
    if not values:
        values = re.findall(r"[^\s]+", block)

    normalized: List[str] = []
    seen: set[str] = set()
    for value in values:
        value = value.strip().strip('"').replace("\\", "/")
        if not value:
            continue
        name = Path(value).name
        if name and name not in seen:
            seen.add(name)
            normalized.append(name)
    return normalized


class LocalBlastDatabaseManager:
    """
    Manage NCBI preformatted BLAST databases and custom FASTA databases.

    The manager never downloads a database implicitly. NCBI's
    update_blastdb.pl is used for official preformatted database downloads,
    while makeblastdb is used for custom FASTA collections.
    """

    def __init__(
        self,
        root: Optional[str | Path] = None,
        executable_manager: Optional[BlastExecutableManager] = None,
    ) -> None:
        configured = os.environ.get(LOCAL_DB_ROOT_ENV, "").strip()
        self.root = Path(root or configured or DEFAULT_LOCAL_DB_ROOT).expanduser()
        self.root = self.root.resolve()
        self.executable_manager = executable_manager or BlastExecutableManager()
        self.manifest_path = self.root / LOCAL_DB_MANIFEST_NAME

    def ensure_root(self) -> Path:
        self.root.mkdir(parents=True, exist_ok=True)
        return self.root

    def _load_manifest(self) -> Dict[str, Any]:
        if not self.manifest_path.is_file():
            return {}
        try:
            payload = json.loads(self.manifest_path.read_text(encoding="utf-8"))
            return payload if isinstance(payload, dict) else {}
        except (OSError, json.JSONDecodeError):
            return {}

    def _save_manifest(self, payload: Dict[str, Any]) -> None:
        self.ensure_root()
        temp_path = self.manifest_path.with_suffix(".tmp")
        temp_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        temp_path.replace(self.manifest_path)

    def _safe_name(self, name: str) -> str:
        """Sanitize a database name while allowing safe subdirectories."""
        raw = str(name).strip().replace("\\", "/")
        if not raw:
            raise BlastAnalysisError("Database name is empty.")

        parts: List[str] = []
        for part in raw.split("/"):
            part = part.strip()
            if not part or part in {".", ".."}:
                if part == "..":
                    raise BlastAnalysisError("Database path traversal is not allowed.")
                continue
            cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", part)
            cleaned = cleaned.strip("._")
            if cleaned in {"", ".", ".."}:
                raise BlastAnalysisError("Invalid database name.")
            parts.append(cleaned)

        if not parts:
            raise BlastAnalysisError("Database name is empty after sanitization.")
        return "/".join(parts)

    def _prefix(self, name: str) -> Path:
        """Return a database prefix guaranteed to remain below self.root."""
        safe_name = self._safe_name(name)
        prefix = (self.root / Path(*safe_name.split("/"))).resolve()
        try:
            prefix.relative_to(self.root)
        except ValueError as exc:
            raise BlastAnalysisError("Database path escapes the configured BLAST root.") from exc
        return prefix

    def _expected_extensions(self, database_type: str) -> Tuple[str, ...]:
        if database_type == "protein":
            return (".pal", ".phr", ".pin", ".psq", ".pdb", ".pog", ".pos", ".pot", ".ptf", ".pto", ".pni")
        return (".nal", ".nhr", ".nin", ".nsq", ".ndb", ".nni", ".nog", ".nos", ".not", ".ntf", ".nto")

    def _find_existing_prefix(self, name: str, database_type: str) -> Path:
        """Resolve a database prefix in flexible local BLAST layouts.

        Supported layouts include:

        1. Flat databases::

               blast_db/testdb.nhr
               blast_db/testdb.nin
               blast_db/testdb.nsq

        2. Database folders::

               blast_db/refseq_rna/refseq_rna.nal
               blast_db/refseq_rna/refseq_rna.00/
               blast_db/refseq_rna/refseq_rna.01/

        3. Nested volume files::

               blast_db/refseq_rna/refseq_rna.00.nhr

        The logical database name remains the requested name; the returned
        path is the physical BLAST prefix used by BLAST+.
        """
        direct = self._prefix(name)
        extensions = self._expected_extensions(database_type)
        if any(Path(str(direct) + extension).is_file() for extension in extensions):
            return direct

        # A database container may have the same name as the logical database,
        # with the actual alias prefix one level below it.
        if direct.is_dir():
            nested_prefix = direct / direct.name
            if any(
                Path(str(nested_prefix) + extension).is_file()
                for extension in extensions
            ):
                return nested_prefix

            # NCBI archives can be extracted volume-by-volume. In that case
            # the first archive may contain the logical alias and shared BLAST
            # metadata inside the ``.00`` volume directory, while ``.01`` ...
            # ``.20`` remain sibling volume directories. The alias is still
            # the logical database; it must not be rejected merely because its
            # physical files live in the first volume directory.
            volume_dir_pattern = re.compile(r"^.+\.\d+$", re.IGNORECASE)
            try:
                alias_candidate = direct / f"{direct.name}.nal"
                protein_alias_candidate = direct / f"{direct.name}.pal"
                if alias_candidate.is_file() or protein_alias_candidate.is_file():
                    return direct

                first_volume = direct / f"{direct.name}.00"
                if first_volume.is_dir():
                    for extension in extensions:
                        candidate = first_volume / f"{direct.name}{extension}"
                        if candidate.is_file():
                            return first_volume / direct.name

                    # More robustly locate the alias in the first numbered
                    # volume directory without treating it as a separate DB.
                    for marker in first_volume.iterdir():
                        if not marker.is_file():
                            continue
                        if marker.name.lower() in {
                            f"{direct.name}.nal".lower(),
                            f"{direct.name}.pal".lower(),
                        }:
                            return first_volume / direct.name
            except OSError:
                pass

        # Search recursively for the requested database prefix. This supports
        # user-chosen roots and nested database containers without hard-coding
        # a Project Synapse path.
        leaf = direct.name
        try:
            matches: List[Path] = []
            volume_dir_pattern = re.compile(r"^.+\.\d+$", re.IGNORECASE)
            for extension in extensions:
                for path in self.root.rglob(f"{leaf}{extension}"):
                    if not path.is_file():
                        continue
                    # Never resolve a logical database prefix from inside a
                    # numbered BLAST volume directory. Those files belong to
                    # the parent database and are not a separate database.
                    if any(volume_dir_pattern.fullmatch(part) for part in path.parent.parts):
                        continue
                    matches.append(path.with_suffix(""))
            if matches:
                matches = sorted(
                    {path.resolve() for path in matches},
                    key=lambda path: (
                        0 if path.parent.name == leaf else 1,
                        len(path.parts),
                        str(path).lower(),
                    ),
                )
                return matches[0]
        except OSError:
            pass

        return direct

    def _volume_container_prefix(self, prefix: Path) -> Path:
        """Return the logical container prefix used to discover sibling volumes.

        Some manually extracted NCBI database archives place the logical alias
        and shared metadata inside ``name.00/`` while the remaining volumes are
        siblings of that directory. In that layout the physical alias prefix is
        ``name.00/name``, but volume discovery must happen relative to
        ``parent/name``.
        """
        base = prefix.name
        parent_name = prefix.parent.name
        if re.fullmatch(rf"{re.escape(base)}\.\d+", parent_name, re.IGNORECASE):
            return prefix.parent.parent / base
        return prefix

    def _discover_volume_prefixes(
        self,
        prefix: Path,
        database_type: str,
    ) -> Dict[str, Path]:
        """Discover BLAST volume prefixes below the alias directory.

        This deliberately handles both ``name.00.nhr`` files and the folder
        layout ``name.00/name.00.nhr`` used by extracted NCBI database archives.
        """
        extensions = self._expected_extensions(database_type)
        base = prefix.name
        volume_pattern = re.compile(rf"^{re.escape(base)}\.(\d+)$", re.IGNORECASE)
        found: Dict[str, Path] = {}

        def add_candidate(candidate: Path) -> None:
            match = volume_pattern.fullmatch(candidate.name)
            if match:
                found[candidate.name] = candidate.resolve()

        container = prefix.parent

        # Flat multi-volume files: name.00.nhr / name.01.nhr / ...
        try:
            for path in container.glob(f"{base}.*"):
                if path.is_file():
                    for extension in extensions:
                        if path.name.lower().endswith(extension.lower()):
                            add_candidate(path.with_suffix(""))
                            break
                elif path.is_dir():
                    add_candidate(path)
        except OSError:
            pass

        # Extracted NCBI archives often create one directory per volume:
        #   refseq_rna.00/refseq_rna.00.nhr
        #   refseq_rna.01/refseq_rna.01.nhr
        try:
            for child in container.iterdir():
                if not child.is_dir() or not volume_pattern.fullmatch(child.name):
                    continue
                volume_prefix = child / child.name
                if any(
                    Path(str(volume_prefix) + extension).is_file()
                    for extension in extensions
                ):
                    add_candidate(volume_prefix)
        except OSError:
            pass

        # Last-resort recursive discovery handles an additional nesting level
        # while remaining restricted to files belonging to this database name.
        try:
            for path in container.rglob("*"):
                if not path.is_file():
                    continue
                if not any(path.name.lower().endswith(ext.lower()) for ext in extensions):
                    continue
                stem = path.name
                for extension in extensions:
                    if stem.lower().endswith(extension.lower()):
                        stem = stem[: -len(extension)]
                        break
                candidate = path.with_name(stem)
                if volume_pattern.fullmatch(candidate.name):
                    add_candidate(candidate)
        except OSError:
            pass

        return dict(sorted(found.items(), key=lambda item: (int(item[0].rsplit('.', 1)[1]), item[0])))

    def _discover_files(
        self,
        name: str,
        database_type: str,
        *,
        prefix: Optional[Path] = None,
    ) -> List[Path]:
        prefix = prefix or self._find_existing_prefix(name, database_type)
        extensions = self._expected_extensions(database_type)
        files: List[Path] = []

        for extension in extensions:
            direct = Path(str(prefix) + extension)
            if direct.is_file():
                files.append(direct)

        volume_discovery_prefix = self._volume_container_prefix(prefix)
        for volume_prefix in self._discover_volume_prefixes(volume_discovery_prefix, database_type).values():
            for extension in extensions:
                candidate = Path(str(volume_prefix) + extension)
                if candidate.is_file():
                    files.append(candidate)

        unique = {path.resolve(): path for path in files}
        return sorted(unique.values(), key=lambda path: str(path).lower())

    def _resolve_expected_volume_prefix(
        self,
        alias_prefix: Path,
        volume_name: str,
        database_type: str,
    ) -> Optional[Path]:
        """Resolve one alias-declared BLAST volume deterministically.

        NCBI extracted archives used by Synapse store each volume in a
        directory such as ``refseq_rna.00/refseq_rna.00.nsq``.  Do not rely
        on recursive discovery for an alias-declared volume: resolve the
        exact volume name first, then accept the flat equivalent as a
        compatibility layout.
        """
        extensions = self._expected_extensions(database_type)
        container = alias_prefix.parent
        volume_dir = container / volume_name
        nested_prefix = volume_dir / volume_name

        if volume_dir.is_dir() and any(
            Path(str(nested_prefix) + extension).is_file()
            for extension in extensions
        ):
            return nested_prefix.resolve()

        flat_prefix = container / volume_name
        if any(
            Path(str(flat_prefix) + extension).is_file()
            for extension in extensions
        ):
            return flat_prefix.resolve()

        # Compatibility with an additional extraction directory level, but
        # only for the exact volume name declared by the alias.
        try:
            for candidate in container.rglob(volume_name):
                if not candidate.is_dir():
                    continue
                nested = candidate / volume_name
                if any(
                    Path(str(nested) + extension).is_file()
                    for extension in extensions
                ):
                    return nested.resolve()
        except OSError:
            pass

        return None


    def _missing_volume_files(
        self,
        volume_prefix: Path,
        database_type: str,
    ) -> List[str]:
        """Return missing mandatory sequence-volume files.

        NCBI multi-volume databases use the sequence core files for each
        volume. Version-5 metadata/index files are not necessarily present
        in every volume, so they are not required for volumes after the first.
        """
        if database_type == "protein":
            required = (".phr", ".pin", ".psq")
        else:
            required = (".nhr", ".nin", ".nsq")
        return [
            extension
            for extension in required
            if not Path(str(volume_prefix) + extension).is_file()
        ]

    def _volume_core_is_complete(
        self,
        volume_prefix: Path,
        database_type: str,
    ) -> bool:
        """Check the mandatory sequence-volume files for one volume."""
        return not self._missing_volume_files(volume_prefix, database_type)


    def scan_database(self, name: str, database_type: Optional[str] = None) -> LocalBlastDatabase:
        safe_name = self._safe_name(name)
        manifest = self._load_manifest()

        metadata = (
            manifest.get(safe_name, {})
            if isinstance(manifest.get(safe_name, {}), dict)
            else {}
        )
        catalog = LOCAL_DATABASE_CATALOG.get(safe_name, {})
        db_type = database_type or metadata.get("database_type") or catalog.get("type")

        if db_type not in {"nucleotide", "protein"}:
            direct_nucl = self._discover_files(safe_name, "nucleotide")
            direct_prot = self._discover_files(safe_name, "protein")
            if direct_prot and not direct_nucl:
                db_type = "protein"
            else:
                db_type = "nucleotide"

        prefix = self._find_existing_prefix(safe_name, db_type)
        files = self._discover_files(safe_name, db_type, prefix=prefix)

        alias_extension = ".pal" if db_type == "protein" else ".nal"
        alias = Path(str(prefix) + alias_extension)

        if db_type == "protein":
            direct_core = all(
                Path(str(prefix) + ext).is_file()
                for ext in (".phr", ".pin", ".psq")
            )
            v5_core = Path(str(prefix) + ".pdb").is_file()
        else:
            direct_core = all(
                Path(str(prefix) + ext).is_file()
                for ext in (".nhr", ".nin", ".nsq")
            )
            v5_core = Path(str(prefix) + ".ndb").is_file()

        volume_discovery_prefix = self._volume_container_prefix(prefix)
        discovered_volume_prefixes = self._discover_volume_prefixes(
            volume_discovery_prefix, db_type
        )
        expected_volume_names = _parse_blast_alias_volume_names(alias)

        # For an alias database, DBLIST is authoritative. Resolve each declared
        # volume by its exact name instead of trusting generic recursive
        # discovery. This is critical for extracted NCBI layouts such as:
        #   refseq_rna/refseq_rna.00/refseq_rna.00.nsq
        #   refseq_rna/refseq_rna.01/refseq_rna.01.nsq
        if expected_volume_names:
            volume_prefixes: Dict[str, Path] = {}
            for volume_name in expected_volume_names:
                resolved_volume = self._resolve_expected_volume_prefix(
                    volume_discovery_prefix, volume_name, db_type
                )

                # IMPORTANT: NCBI BLAST aliases may use an unpadded numeric
                # volume token (``.0``, ``.1`` ...) while the extracted archive
                # is stored with two digits (``.00``, ``.01`` ...). Treat both
                # representations as the same volume. Never construct a new
                # filesystem name by dropping zero padding. Resolve against
                # the actually discovered directory/file prefix instead.
                if resolved_volume is None:
                    expected_index = _extract_blast_volume_index(volume_name)
                    if expected_index is not None:
                        for observed_name, observed_prefix in discovered_volume_prefixes.items():
                            observed_index = _extract_blast_volume_index(observed_name)
                            if observed_index == expected_index:
                                resolved_volume = observed_prefix
                                break

                if resolved_volume is not None:
                    volume_prefixes[volume_name] = resolved_volume
        else:
            volume_prefixes = discovered_volume_prefixes

        observed_volume_names = list(volume_prefixes)

        missing_volumes = [
            volume
            for volume in expected_volume_names
            if volume not in volume_prefixes
        ]
        missing_volume_files: Dict[str, List[str]] = {}
        incomplete_volumes: List[str] = []
        for volume, volume_prefix in volume_prefixes.items():
            missing = self._missing_volume_files(volume_prefix, db_type)
            if missing:
                incomplete_volumes.append(volume)
                missing_volume_files[volume] = missing

        if expected_volume_names:
            all_expected_present = not missing_volumes
            all_expected_complete = all(
                volume in volume_prefixes
                and self._volume_core_is_complete(volume_prefixes[volume], db_type)
                for volume in expected_volume_names
            )
            structural = alias.is_file() and all_expected_present and all_expected_complete
        elif volume_prefixes:
            structural = alias.is_file() and not incomplete_volumes
        else:
            structural = alias.is_file() or direct_core or v5_core

        modified = None
        size = 0
        existing_files = [path for path in files if path.exists()]
        if existing_files:
            size = sum(path.stat().st_size for path in existing_files)
            newest = max(path.stat().st_mtime for path in existing_files)
            modified = datetime.fromtimestamp(newest, timezone.utc).replace(
                microsecond=0
            ).isoformat()

        # A catalog database such as refseq_rna may physically live in
        # blast_db/refseq_rna/refseq_rna. Keep the logical name stable while
        # retaining the physical prefix in ``path``. For explicitly nested
        # custom names (e.g. local_test/synapse_test), preserve the nested name.
        resolved_name = safe_name
        try:
            resolved_relative = prefix.relative_to(self.root)
            if "/" in safe_name:
                resolved_name = resolved_relative.as_posix()
            else:
                parent_name = prefix.parent.name
                is_numbered_volume_parent = bool(
                    re.fullmatch(rf"{re.escape(safe_name)}\.\d+", parent_name, re.IGNORECASE)
                )
                is_standard_catalog_prefix = (
                    prefix.parent.name == safe_name and prefix.name == safe_name
                )
                if not is_numbered_volume_parent and not is_standard_catalog_prefix:
                    if str(resolved_relative.parent) != ".":
                        resolved_name = resolved_relative.as_posix()
        except ValueError:
            pass

        resolved_metadata = manifest.get(resolved_name, {})
        if isinstance(resolved_metadata, dict) and resolved_metadata:
            metadata = resolved_metadata

        if expected_volume_names:
            volume_status = (
                f"{len(expected_volume_names) - len(missing_volumes)}/"
                f"{len(expected_volume_names)} volumes present"
            )
        elif volume_prefixes:
            volume_status = f"{len(volume_prefixes)} volume(s) detected"
        else:
            volume_status = "Single-volume database"

        if structural:
            verification_status = "STRUCTURALLY_VALID"
            verification_message = (
                f"BLAST database structure detected. {volume_status}. "
                "Use blastdbcmd for definitive verification."
            )
        elif existing_files:
            verification_status = "INCOMPLETE"
            details: List[str] = []
            if missing_volumes:
                details.append("missing: " + ", ".join(missing_volumes))
            if incomplete_volumes:
                details.append("incomplete: " + ", ".join(incomplete_volumes))
            verification_message = (
                "BLAST database files were found, but the database is incomplete. "
                + ("; ".join(details) if details else volume_status)
            )
        else:
            verification_status = "NOT_INSTALLED"
            verification_message = "Required BLAST database files were not detected."

        return LocalBlastDatabase(
            name=resolved_name,
            database_type=db_type,
            path=str(prefix),
            source=metadata.get("source") or catalog.get("source") or "Local",
            installed=bool(existing_files),
            structurally_valid=structural,
            verification_status=verification_status,
            verification_message=verification_message,
            file_count=len(existing_files),
            size_bytes=size,
            modified_at=modified,
            volumes=observed_volume_names,
            volume_paths=[
                str(path.parent if path.parent.name == path.name else path)
                for path in volume_prefixes.values()
            ],
            volume_count=len(volume_prefixes),
            expected_volume_count=(len(expected_volume_names) if expected_volume_names else None),
            missing_volumes=missing_volumes,
            incomplete_volumes=incomplete_volumes,
            missing_volume_files=missing_volume_files,
            manifest_metadata=metadata,
        )


    def list_databases(self) -> List[LocalBlastDatabase]:
        manifest = self._load_manifest()
        names = set(LOCAL_DATABASE_CATALOG)
        names.update(
            name for name in manifest
            if name != LOCAL_DB_MANIFEST_NAME and isinstance(name, str)
        )

        # Discover custom/nested databases, but NEVER promote a BLAST volume
        # directory (for example refseq_rna.00/) to a database of its own.
        # NCBI multi-volume databases are represented by one logical alias
        # (refseq_rna.nal) plus refseq_rna.00, refseq_rna.01, ... volumes.
        try:
            patterns = ("*.nal", "*.pal", "*.ndb", "*.pdb")
            volume_dir_pattern = re.compile(r"^.+\.\d+$", re.IGNORECASE)

            for pattern in patterns:
                for marker in self.root.rglob(pattern):
                    if not marker.is_file():
                        continue
                    try:
                        relative = marker.relative_to(self.root)
                    except ValueError:
                        continue

                    parts = relative.parts
                    # Ignore any marker physically stored inside a numbered
                    # volume directory. Such files belong to the parent
                    # database and must not create a second database entry.
                    if any(volume_dir_pattern.fullmatch(part) for part in parts[:-1]):
                        continue

                    stem = marker.stem
                    parent = marker.parent

                    # A volume file such as refseq_rna.00.ndb must also never
                    # become a standalone database.
                    if re.fullmatch(r".+\.\d+", stem):
                        continue

                    # Accept the supported container layout:
                    #   blast_db/name/name.nal
                    # as the logical database name ``name`` rather than the
                    # physical prefix ``name/name``. Flat databases and
                    # genuinely nested custom databases keep their relative path.
                    if parent.name == stem and parent != self.root:
                        try:
                            name = parent.relative_to(self.root).as_posix()
                        except ValueError:
                            name = stem
                    else:
                        name = relative.with_suffix("").as_posix()

                    if name:
                        names.add(name)
        except OSError:
            pass

        databases = [self.scan_database(name) for name in sorted(names, key=str.lower)]

        # Deduplicate aliases that resolve to the same physical prefix.
        # Prefer the shortest logical/catalog name (for example ``refseq_rna``)
        # over accidental nested aliases such as ``refseq_rna/refseq_rna.00/...``.
        unique: Dict[str, LocalBlastDatabase] = {}
        for database in databases:
            key = str(Path(database.path).resolve()).lower()
            current = unique.get(key)
            if current is None:
                unique[key] = database
                continue
            current_parts = len(Path(current.name).parts)
            candidate_parts = len(Path(database.name).parts)
            if candidate_parts < current_parts:
                unique[key] = database
            elif candidate_parts == current_parts and database.name.lower() < current.name.lower():
                unique[key] = database
        return sorted(unique.values(), key=lambda item: item.name.lower())


    def database_exists(self, name: str) -> bool:
        return self.scan_database(name).exists

    def _runtime_alias_parent(self, database: LocalBlastDatabase) -> Path:
        """Return a stable directory in which a runtime alias can resolve volumes.

        For nested NCBI layouts the physical volume prefixes are siblings below
        the database container, so the temporary alias must live in that same
        container (or a child of it). This keeps DBLIST entries relative and
        independent of Windows paths containing spaces.
        """
        database_path = Path(database.path).expanduser().resolve()
        volume_directory = database_path.parent
        if re.fullmatch(r".+\.\d+", volume_directory.name):
            return volume_directory.parent
        return volume_directory

    def verify_database(self, name: str, run_blastdbcmd: bool = True) -> LocalBlastDatabase:
        database = self.scan_database(name)
        if not database.structurally_valid:
            database.verification_status = "FAILED"
            details: List[str] = []
            if database.missing_volumes:
                details.append("missing volumes: " + ", ".join(database.missing_volumes))
            if database.missing_volume_files:
                file_details = "; ".join(
                    f"{volume} [{', '.join(files)}]"
                    for volume, files in database.missing_volume_files.items()
                )
                details.append("missing files: " + file_details)
            database.verification_message = (
                "Database structure is incomplete. "
                + ("; ".join(details) if details else "No complete database core was detected.")
            )
            return database

        if not run_blastdbcmd:
            database.verification_status = "STRUCTURALLY_VALID"
            database.verification_message = "File structure is valid."
            return database

        try:
            blastdbcmd = self.executable_manager.require("blastdbcmd")
        except BlastAnalysisError as exc:
            database.verification_status = "STRUCTURALLY_VALID_TOOL_MISSING"
            database.verification_message = (
                "File structure is valid, but blastdbcmd is unavailable. "
                f"{exc}"
            )
            return database

        runtime_directory: Optional[tempfile.TemporaryDirectory[str]] = None
        try:
            if _needs_runtime_blast_alias(database):
                runtime_parent = self._runtime_alias_parent(database)
                runtime_directory = tempfile.TemporaryDirectory(
                    prefix="synapse_blast_verify_",
                    dir=str(runtime_parent),
                )
                database_reference, database_working_directory = _create_runtime_blast_alias(
                    database, runtime_directory.name
                )
            else:
                database_reference, database_working_directory = (
                    _build_local_database_command_reference(database.path)
                )

            if _needs_runtime_blast_alias(database):
                env = _build_local_blast_environment(
                    self._runtime_alias_parent(database)
                )
            else:
                env = _build_local_blast_environment(database.path)
            command = [blastdbcmd, "-db", database_reference, "-info"]
            completed = subprocess.run(
                command,
                cwd=str(database_working_directory),
                env=env,
                capture_output=True,
                text=True,
                timeout=120,
                check=False,
            )
        except (OSError, subprocess.SubprocessError, BlastAnalysisError) as exc:
            database.verification_status = "VERIFICATION_ERROR"
            database.verification_message = str(exc)
            return database
        finally:
            if runtime_directory is not None:
                runtime_directory.cleanup()

        output = (completed.stdout or "") + "\r\n" + (completed.stderr or "")
        if completed.returncode == 0:
            parsed = _parse_blastdbcmd_info(output)
            database.title = parsed["title"]
            database.sequence_count = parsed["sequence_count"]
            database.total_length = parsed["total_length"]
            database.blastdb_version = parsed["blastdb_version"]
            database.volumes = parsed["volumes"]
            database.verification_status = "VERIFIED"
            database.verification_message = output.strip() or "blastdbcmd verification succeeded."
        else:
            database.verification_status = "FAILED"
            database.verification_message = output.strip() or (
                f"blastdbcmd exited with code {completed.returncode}."
            )
        return database

    def update_ncbi_database(
        self,
        name: str,
        *,
        timeout: Optional[float] = None,
    ) -> LocalBlastDatabase:
        if name not in LOCAL_DATABASE_CATALOG:
            raise BlastAnalysisError(
                f"'{name}' is not in the Synapse NCBI database catalog."
            )

        script = self.executable_manager.resolve_update_script()
        if not script:
            raise BlastAnalysisError(
                "update_blastdb.pl was not found. Install the NCBI BLAST+ package "
                "and make sure its bin directory is discoverable."
            )

        perl = self.executable_manager.resolve_perl()
        if not perl:
            raise BlastAnalysisError(
                "Perl was not found. NCBI's update_blastdb.pl requires Perl."
            )

        self.ensure_root()
        command = [
            perl,
            script,
            "--passive",
            "--decompress",
            name,
        ]

        try:
            completed = subprocess.run(
                command,
                cwd=str(self.root),
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise BlastAnalysisError(
                f"Database update timed out after {timeout} seconds."
            ) from exc
        except OSError as exc:
            raise BlastAnalysisError(f"Could not start database update: {exc}") from exc

        output = ((completed.stdout or "") + "\r\n" + (completed.stderr or "")).strip()
        if completed.returncode != 0:
            raise BlastAnalysisError(
                f"NCBI database update failed with exit code {completed.returncode}.\r\n"
                f"{output[-4000:]}"
            )

        manifest = self._load_manifest()
        manifest[name] = {
            "database_type": LOCAL_DATABASE_CATALOG[name]["type"],
            "source": LOCAL_DATABASE_CATALOG[name]["source"],
            "description": LOCAL_DATABASE_CATALOG[name]["description"],
            "last_update_command": "update_blastdb.pl --passive --decompress " + name,
            "updated_at_utc": _utc_now(),
        }
        self._save_manifest(manifest)

        database = self.scan_database(name)
        database.verification_message = output[-4000:] or "Database update completed."
        return database

    def create_custom_database(
        self,
        fasta_path: str | Path,
        name: str,
        database_type: str,
        *,
        title: Optional[str] = None,
        overwrite: bool = False,
    ) -> LocalBlastDatabase:
        database_type = database_type.lower().strip()
        if database_type not in {"nucleotide", "protein"}:
            raise BlastAnalysisError("database_type must be 'nucleotide' or 'protein'.")

        fasta = Path(fasta_path).expanduser().resolve()
        if not fasta.is_file():
            raise BlastAnalysisError(f"FASTA file not found: {fasta}")

        safe_name = self._safe_name(name)
        existing = self.scan_database(safe_name, database_type)
        if existing.installed and not overwrite:
            raise BlastAnalysisError(
                f"Database '{safe_name}' already exists. "
                "Choose another name or enable overwrite."
            )

        makeblastdb = self.executable_manager.require("makeblastdb")
        self.ensure_root()
        prefix = self._prefix(safe_name)
        dbtype_flag = "prot" if database_type == "protein" else "nucl"

        command = [
            makeblastdb,
            "-in",
            str(fasta),
            "-dbtype",
            dbtype_flag,
            "-parse_seqids",
            "-out",
            str(prefix),
            "-title",
            title or safe_name,
        ]

        try:
            completed = subprocess.run(
                command,
                cwd=str(self.root),
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise BlastAnalysisError(f"Could not start makeblastdb: {exc}") from exc

        output = ((completed.stdout or "") + "\r\n" + (completed.stderr or "")).strip()
        if completed.returncode != 0:
            raise BlastAnalysisError(
                f"makeblastdb failed with exit code {completed.returncode}.\r\n"
                f"{output[-4000:]}"
            )

        manifest = self._load_manifest()
        manifest[safe_name] = {
            "database_type": database_type,
            "source": "Custom FASTA",
            "description": title or f"Custom FASTA database: {fasta.name}",
            "fasta_source": str(fasta),
            "created_at_utc": _utc_now(),
        }
        self._save_manifest(manifest)

        database = self.scan_database(safe_name, database_type)
        database.verification_message = output[-4000:] or "Custom database created."
        return database

    def remove_database(self, name: str) -> List[str]:
        safe_name = self._safe_name(name)
        database = self.scan_database(safe_name)
        removed: List[str] = []

        if not database.installed:
            return removed

        prefix = Path(database.path).resolve()
        database_type = database.database_type
        files = self._discover_files(safe_name, database_type, prefix=prefix)
        volume_prefixes = self._discover_volume_prefixes(prefix, database_type)

        # Remove only artifacts belonging to this database. Archive files,
        # unrelated databases and the parent database container are preserved.
        targets: set[Path] = {path.resolve() for path in files if path.is_file()}
        targets.update(path.resolve() for path in volume_prefixes.values() if path.is_dir())

        # If the logical database prefix itself is a directory (unusual but
        # supported), remove only its database-owned files, never the parent.
        for target in sorted(targets, key=lambda path: (len(path.parts), str(path).lower()), reverse=True):
            try:
                if target.is_file():
                    target.unlink()
                    removed.append(str(target))
                elif target.is_dir():
                    target.rmdir()
                    removed.append(str(target))
            except OSError as exc:
                raise BlastAnalysisError(f"Could not remove '{target}': {exc}") from exc

        manifest = self._load_manifest()
        for manifest_name in {safe_name, database.name}:
            if manifest_name in manifest:
                del manifest[manifest_name]
        self._save_manifest(manifest)

        return removed


    def database_info(self, name: str) -> Dict[str, Any]:
        database = self.scan_database(name)
        return {
            "name": database.name,
            "database_type": database.database_type,
            "path": database.path,
            "source": database.source,
            "installed": database.installed,
            "structurally_valid": database.structurally_valid,
            "verification_status": database.verification_status,
            "verification_message": database.verification_message,
            "file_count": database.file_count,
            "size_bytes": database.size_bytes,
            "size_mb": round(database.size_mb, 2),
            "size_gb": round(database.size_gb, 3),
            "modified_at": database.modified_at,
            "title": database.title,
            "sequence_count": database.sequence_count,
            "total_length": database.total_length,
            "blastdb_version": database.blastdb_version,
            "volumes": database.volumes,
            "volume_paths": database.volume_paths,
            "volume_count": database.volume_count,
            "expected_volume_count": database.expected_volume_count,
            "missing_volumes": database.missing_volumes,
            "incomplete_volumes": database.incomplete_volumes,
            "missing_volume_files": database.missing_volume_files,
            "manifest_metadata": database.manifest_metadata,
        }


def _build_local_database_command_reference(database_path: str | Path) -> Tuple[str, Path]:
    """Return a BLAST database reference and working directory.

    BLAST+ on Windows is generally path-safe when launched through
    ``subprocess`` with an argument list, but some BLAST+ builds and database
    layouts behave more reliably when the database prefix is supplied relative
    to its own directory. This also avoids ambiguity for paths containing
    spaces (for example ``Project Synapse``).

    The returned reference is only the database prefix name; the caller must
    set ``cwd`` to the returned parent directory. No shell is used.
    """
    resolved = Path(database_path).expanduser().resolve()
    parent = resolved.parent
    return resolved.name, parent


def _volume_prefix_from_stored_path(volume_path: str | Path) -> Path:
    """Resolve a stored volume path to the actual BLAST database prefix.

    Nested NCBI extraction layouts are stored by Synapse as directories for
    display purposes, e.g. ``.../refseq_rna.01``. The BLAST prefix inside that
    directory is ``.../refseq_rna.01/refseq_rna.01``. Flat prefixes are kept as
    they are.
    """
    resolved = Path(volume_path).expanduser().resolve()
    if resolved.is_dir():
        nested_prefix = resolved / resolved.name
        extensions = (".nhr", ".nin", ".nsq", ".phr", ".pin", ".psq")
        if any(Path(str(nested_prefix) + extension).is_file() for extension in extensions):
            return nested_prefix
    return resolved


def _create_runtime_blast_alias(
    database: LocalBlastDatabase,
    runtime_directory: str | Path,
) -> Tuple[str, Path]:
    """Create a temporary flat BLAST view for a nested multi-volume database.

    NCBI's preformatted multi-volume databases are normally exposed to BLAST
    as a flat set of files:

        refseq_rna.nal
        refseq_rna.00.nhr / .nin / .nsq / ...
        refseq_rna.01.nhr / .nin / .nsq / ...
        ...

    Synapse also supports the extraction layout produced by manually unpacking
    each NCBI archive into a volume directory, for example::

        refseq_rna/refseq_rna.00/refseq_rna.00.nhr
        refseq_rna/refseq_rna.01/refseq_rna.01.nhr
        ...

    BLAST+ does not treat that nested layout as the canonical multi-volume
    database when it is presented through a synthetic alias. In particular,
    v5 metadata such as ``refseq_rna.ndb`` is then interpreted incorrectly for
    later volumes. The robust solution is therefore a temporary *flat view*.

    The view uses hard links, not file copies. Hard links consume essentially
    no additional database storage and leave the original NCBI files untouched.
    The runtime directory is on the same filesystem as the source database so
    Windows can create the hard links. The alias keeps the canonical database
    name and volume names, which lets BLAST+ use its normal multi-volume logic.
    """
    requested_runtime = Path(runtime_directory).expanduser().resolve()
    database_path = Path(database.path).expanduser().resolve()
    requested_runtime.mkdir(parents=True, exist_ok=True)

    extension = ".pal" if database.database_type == "protein" else ".nal"
    logical_name = Path(database.name).name
    logical_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", logical_name) or "synapse_database"

    # A BLAST multi-volume alias must use the same logical root as its volume
    # prefixes (for example refseq_rna + refseq_rna.00 ... refseq_rna.20).
    # This is intentionally NOT a *_runtime.nal alias.
    alias_path = requested_runtime / f"{logical_name}{extension}"

    volume_prefixes = [
        _volume_prefix_from_stored_path(path)
        for path in database.volume_paths
    ]
    if not volume_prefixes:
        raise BlastAnalysisError(
            f"Cannot create runtime BLAST view for '{database.name}': "
            "no volume paths were detected."
        )

    # The source database must be on the same filesystem because Windows
    # hard-links cannot cross volumes. tempfile.TemporaryDirectory callers
    # already place the directory next to the database for this reason.
    try:
        requested_runtime.relative_to(database_path.parent.parent)
    except ValueError:
        # This is only a diagnostic guard; a same-volume check is performed by
        # os.link below and produces the actionable error if needed.
        pass

    db_entries: List[str] = []
    created_links: List[Path] = []

    try:
        for volume_prefix in volume_prefixes:
            volume_name = volume_prefix.name
            if not re.fullmatch(rf"{re.escape(logical_name)}\.\d+", volume_name, re.IGNORECASE):
                raise BlastAnalysisError(
                    f"Unexpected BLAST volume name '{volume_name}' for database '{logical_name}'."
                )

            # Reproduce the canonical flat BLAST volume prefix. Link every
            # source file belonging to this volume, including v5 metadata on
            # the first volume when present.
            source_parent = volume_prefix.parent
            pattern = f"{volume_name}.*"
            source_files = sorted(source_parent.glob(pattern))
            if not source_files:
                raise BlastAnalysisError(
                    f"No files found for BLAST volume '{volume_prefix}'."
                )

            for source_file in source_files:
                if not source_file.is_file():
                    continue
                target = requested_runtime / source_file.name
                try:
                    if target.exists() or target.is_symlink():
                        target.unlink()
                    os.link(str(source_file), str(target))
                    created_links.append(target)
                except OSError as exc:
                    raise BlastAnalysisError(
                        "Could not create a temporary hard link for BLAST database "
                        f"file '{source_file}': {exc}"
                    ) from exc

            db_entries.append(volume_name)

        # The original NCBI alias contains the database-level v5 files in the
        # first volume directory. Those files do not carry a volume suffix and
        # therefore are not covered by the volume glob above.
        first_volume_dir = volume_prefixes[0].parent
        for source_file in sorted(first_volume_dir.iterdir()):
            if not source_file.is_file():
                continue
            if source_file.name == alias_path.name:
                continue
            # Only link database-level metadata that BLAST may need. Taxonomy
            # files are included as well so the runtime view remains complete.
            if re.fullmatch(rf"{re.escape(logical_name)}\.\d+\..+", source_file.name, re.IGNORECASE):
                continue
            if source_file.name.lower() in {
                "taxdb.btd",
                "taxdb.bti",
                "taxonomy4blast.sqlite3",
            } or source_file.name.lower() in {
                f"{logical_name}.ndb",
                f"{logical_name}.nos",
                f"{logical_name}.not",
                f"{logical_name}.ntf",
                f"{logical_name}.nto",
            }:
                target = requested_runtime / source_file.name
                try:
                    if target.exists() or target.is_symlink():
                        target.unlink()
                    os.link(str(source_file), str(target))
                    created_links.append(target)
                except OSError as exc:
                    raise BlastAnalysisError(
                        f"Could not create a temporary hard link for '{source_file}': {exc}"
                    ) from exc

        title = database.title or f"Synapse runtime view for {database.name}"
        alias_path.write_text(
            "# Synapse temporary BLAST alias\r\n"
            f"TITLE {title}\r\n"
            f"DBLIST {' '.join(db_entries)}\r\n",
            encoding="utf-8",
            newline="",
        )

        return alias_path.stem, alias_path.parent

    except Exception:
        # Never leave a partially constructed runtime view behind.
        for path in reversed(created_links):
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass
        try:
            alias_path.unlink(missing_ok=True)
        except OSError:
            pass
        raise

def _needs_runtime_blast_alias(database: LocalBlastDatabase) -> bool:
    """Return whether the physical database layout needs a runtime alias."""
    if database.volume_count <= 1 or not database.volume_paths:
        return False

    database_path = Path(database.path).expanduser().resolve()
    parent_name = database_path.parent.name
    return bool(re.fullmatch(r".+\.\d+", parent_name))


def _build_local_blast_environment(database_path: str | Path) -> Dict[str, str]:
    """Build a BLAST+ environment with an explicit database search path."""
    env = os.environ.copy()
    resolved = Path(database_path).expanduser().resolve()
    parent = resolved if resolved.is_dir() else resolved.parent
    existing = env.get("BLASTDB", "").strip()
    parent_text = str(parent)
    if existing:
        env["BLASTDB"] = os.pathsep.join([parent_text, existing])
    else:
        env["BLASTDB"] = parent_text
    return env


def _should_run_short_query_rescue(
    query: BlastQuery,
    parameters: BlastParameters,
) -> bool:
    """Return whether Synapse may add one automatic short-query rescue run."""
    return (
        parameters.program.lower().strip() == "blastn"
        and query.sequence_type == "DNA"
        and query.length <= 200
        and parameters.word_size is None
        and parameters.task is None
        and parameters.filter_query is None
    )


class LocalBlastEngine:
    """Execute BLAST+ locally against an installed database."""

    def __init__(
        self,
        executable_manager: Optional[BlastExecutableManager] = None,
    ) -> None:
        self.executable_manager = executable_manager or BlastExecutableManager()

    def run(
        self,
        query: BlastQuery,
        parameters: BlastParameters,
        database: LocalBlastDatabase,
        *,
        timeout: float = DEFAULT_MAX_WAIT_SECONDS,
        allow_sensitive_retry: bool = True,
    ) -> BlastRun:
        if not database.exists:
            raise BlastAnalysisError(
                f"Local BLAST database '{database.name}' is not installed or incomplete."
            )

        expected_type = PROGRAMS.get(parameters.program.lower().strip(), {}).get(
            "database_type"
        )
        if expected_type and expected_type != database.database_type:
            raise BlastAnalysisError(
                f"{parameters.program} requires a {expected_type} database, "
                f"but '{database.name}' is {database.database_type}."
            )

        program = parameters.program.lower().strip()
        if program not in PROGRAMS:
            raise BlastAnalysisError(f"Unsupported BLAST program: {program}")

        blast_executable = self.executable_manager.require(program)

        # Validate query type against the selected BLAST program.
        expected_query = PROGRAMS[program]["query_type"]
        if expected_query != query.sequence_type:
            raise BlastAnalysisError(
                f"{program} expects a {expected_query} query, "
                f"but the query is {query.sequence_type}."
            )

        local_parameters = BlastParameters(
            program=program,
            database=database.name,
            evalue=parameters.evalue,
            max_target_seqs=parameters.max_target_seqs,
            word_size=parameters.word_size,
            gap_open=parameters.gap_open,
            gap_extend=parameters.gap_extend,
            matrix=parameters.matrix,
            filter_query=parameters.filter_query,
            num_threads=parameters.num_threads,
            task=parameters.task,
            tool=parameters.tool,
            email=parameters.email,
        )

        with tempfile.TemporaryDirectory(prefix="synapse_blast_") as temp_dir:
            temp = Path(temp_dir)
            query_file = temp / "query.fasta"
            output_file = temp / "result.xml"

            query_file.write_text(
                f">{query.name}\r\n{query.sequence}\r\n",
                encoding="utf-8",
            )

            if _needs_runtime_blast_alias(database):
                runtime_parent = Path(database.path).expanduser().resolve().parent
                if re.fullmatch(r".+\.\d+", runtime_parent.name):
                    runtime_parent = runtime_parent.parent
                runtime_directory = Path(
                    tempfile.mkdtemp(
                        prefix="synapse_blast_runtime_",
                        dir=str(runtime_parent),
                    )
                )
                database_reference, database_working_directory = _create_runtime_blast_alias(
                    database, runtime_directory
                )
            else:
                database_reference, database_working_directory = (
                    _build_local_database_command_reference(database.path)
                )

            command = [
                blast_executable,
                "-query",
                str(query_file),
                "-db",
                database_reference,
                "-out",
                str(output_file),
                "-outfmt",
                "5",
                "-evalue",
                str(local_parameters.evalue),
                "-max_target_seqs",
                str(local_parameters.max_target_seqs),
            ]

            if local_parameters.num_threads is not None:
                if local_parameters.num_threads < 1:
                    raise BlastAnalysisError("num_threads must be at least 1.")
                command.extend(["-num_threads", str(local_parameters.num_threads)])

            if local_parameters.task:
                command.extend(["-task", local_parameters.task])

            if local_parameters.word_size is not None:
                command.extend(["-word_size", str(local_parameters.word_size)])
            if local_parameters.gap_open is not None:
                command.extend(["-gapopen", str(local_parameters.gap_open)])
            if local_parameters.gap_extend is not None:
                command.extend(["-gapextend", str(local_parameters.gap_extend)])
            if local_parameters.matrix:
                command.extend(["-matrix", local_parameters.matrix])
            if local_parameters.filter_query:
                command.extend(["-seg" if program in {"blastp", "blastx", "tblastn", "tblastx"} else "-dust",
                                local_parameters.filter_query])

            if _needs_runtime_blast_alias(database):
                env = _build_local_blast_environment(
                    runtime_parent
                )
            else:
                env = _build_local_blast_environment(database.path)

            started = _utc_now()
            try:
                completed = subprocess.run(
                    command,
                    cwd=str(database_working_directory),
                    env=env,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    check=False,
                )
            except subprocess.TimeoutExpired as exc:
                raise BlastAnalysisError(
                    f"Local BLAST search exceeded the timeout of {timeout:.0f} seconds."
                ) from exc
            except OSError as exc:
                raise BlastAnalysisError(f"Could not start local BLAST: {exc}") from exc

            stderr = (completed.stderr or "").strip()
            stdout = (completed.stdout or "").strip()
            if completed.returncode != 0:
                detail = stderr or stdout or "No diagnostic output was returned."
                raise BlastAnalysisError(
                    f"Local {program} failed with exit code {completed.returncode}.\r\n"
                    f"{detail[-4000:]}"
                )

            if not output_file.is_file():
                raise BlastAnalysisError(
                    "Local BLAST completed without producing an XML result file."
                )

            xml_text = output_file.read_text(encoding="utf-8", errors="replace")

        rid = f"LOCAL-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
        completed_at = _utc_now()
        run = parse_blast_xml(
            xml_text,
            query=query,
            rid=rid,
            parameters=local_parameters,
            submitted_at=started,
            completed_at=completed_at,
            database_type_override=database.database_type,
        )
        run.metadata.engine = "Local BLAST+"
        run.metadata.local_database_path = database.path
        run.metadata.ncbi_status = "LOCAL_COMPLETE"
        if stderr:
            run.metadata.warnings.append(stderr[-2000:])
        run.interpretation_notes.append(
            "This result was generated by the locally installed NCBI BLAST+ executable."
        )
        run.interpretation_notes.append(
            "The local database is a specific installed snapshot; database contents "
            "may differ from the current NCBI remote database."
        )

        # Short DNA queries can miss local similarities with the standard BLASTN
        # seed model. For an otherwise default local BLASTN search, Synapse performs
        # one conservative sensitivity rescue using NCBI's blastn-short task. The
        # reporting E-value is kept unchanged, so the rescue does not silently turn
        # weak matches into significant hits. Explicit user parameters are never
        # overridden.
        short_query_eligible = _should_run_short_query_rescue(query, parameters)

        if allow_sensitive_retry and short_query_eligible and not run.hits:
            rescue_parameters = BlastParameters(
                program=program,
                database=database.name,
                evalue=parameters.evalue,
                max_target_seqs=parameters.max_target_seqs,
                word_size=None,
                gap_open=parameters.gap_open,
                gap_extend=parameters.gap_extend,
                matrix=parameters.matrix,
                filter_query="no",
                num_threads=parameters.num_threads,
                task="blastn-short",
                tool=parameters.tool,
                email=parameters.email,
            )

            rescue_run = self.run(
                query,
                rescue_parameters,
                database,
                timeout=timeout,
                allow_sensitive_retry=False,
            )
            rescue_run.interpretation_notes.insert(
                0,
                "Standard BLASTN returned no hits; Synapse automatically performed "
                "one short-query sensitivity rescue using the blastn-short task.",
            )
            rescue_run.interpretation_notes.insert(
                1,
                "The sensitivity rescue retained the original E-value threshold and "
                "did not override explicit user search parameters.",
            )
            if rescue_run.hits:
                rescue_run.interpretation_notes.insert(
                    2,
                    "Displayed hits are from the sensitive short-query search and "
                    "should be interpreted using identity, coverage and E-value together.",
                )
            else:
                rescue_run.interpretation_notes.insert(
                    2,
                    "No hits were found in either the standard or short-query search.",
                )
            return rescue_run

        return run


# ============================================================================
# BLAST XML PARSER
# ============================================================================

def _parse_hsp(hsp_element: ET.Element) -> BlastHSP:
    return BlastHSP(
        number=_to_int(_child_text(hsp_element, "Hsp_num"), 1) or 1,
        bit_score=_to_float(_child_text(hsp_element, "Hsp_bit-score"), 0.0) or 0.0,
        raw_score=_to_float(_child_text(hsp_element, "Hsp_score")),
        evalue=_to_float(_child_text(hsp_element, "Hsp_evalue"), math.inf) or math.inf,
        query_start=_to_int(_child_text(hsp_element, "Hsp_query-from"), 0) or 0,
        query_end=_to_int(_child_text(hsp_element, "Hsp_query-to"), 0) or 0,
        subject_start=_to_int(_child_text(hsp_element, "Hsp_hit-from"), 0) or 0,
        subject_end=_to_int(_child_text(hsp_element, "Hsp_hit-to"), 0) or 0,
        query_frame=_to_int(_child_text(hsp_element, "Hsp_query-frame")),
        subject_frame=_to_int(_child_text(hsp_element, "Hsp_hit-frame")),
        identity=_to_int(_child_text(hsp_element, "Hsp_identity"), 0) or 0,
        positive=_to_int(_child_text(hsp_element, "Hsp_positive"), 0) or 0,
        gaps=_to_int(_child_text(hsp_element, "Hsp_gaps"), 0) or 0,
        alignment_length=_to_int(_child_text(hsp_element, "Hsp_align-len"), 0) or 0,
        query_sequence=_child_text(hsp_element, "Hsp_qseq"),
        subject_sequence=_child_text(hsp_element, "Hsp_hseq"),
        midline=_child_text(hsp_element, "Hsp_midline"),
    )


def _infer_scientific_name_from_definition(definition: str) -> Optional[str]:
    """Infer a binomial organism name when BLAST XML omits taxonomy fields.

    Some NCBI BLAST XML responses provide a useful ``Hit_def`` but no
    ``Hit_sciname``/``Hit_taxid`` element. In that case Synapse may extract a
    conservative binomial-looking name from the definition for display. This
    is explicitly an inference from the returned definition, not a taxonomy
    lookup and not an authoritative NCBI taxonomy assignment.
    """
    text = re.sub(r"\s+", " ", (definition or "").strip())
    if not text:
        return None

    bracketed = re.findall(
        r"\[\s*([A-Z][A-Za-z-]{2,})\s+([a-z][A-Za-z-]{2,})\s*\]",
        text,
    )
    if bracketed:
        genus, species = bracketed[0]
        return f"{genus} {species}"

    candidates = re.findall(
        r"(?<![A-Za-z-])([A-Z][A-Za-z-]{2,})\s+([a-z][A-Za-z-]{2,})(?![A-Za-z-])",
        text,
    )
    if not candidates:
        return None

    excluded = {
        "Predicted",
        "Putative",
        "Hypothetical",
        "Uncharacterized",
        "Unnamed",
        "Novel",
        "Example",
        "Reference",
        "Sequence",
        "Protein",
        "Transcript",
    }
    for genus, species in candidates:
        if genus not in excluded:
            return f"{genus} {species}"

    return None


def _taxonomy_source(hit: BlastHit) -> str:
    """Return how Synapse obtained the displayed taxonomy label."""
    if hit.taxonomy_ids:
        return "NCBI BLAST taxonomy fields"
    if hit.scientific_names:
        return "NCBI BLAST scientific-name field"
    if _infer_scientific_name_from_definition(hit.definition):
        return "inferred from BLAST hit definition"
    return "not available"


def _parse_hit(hit_element: ET.Element) -> BlastHit:
    hsp_container = next(
        (node for node in hit_element if _local_tag(node) == "Hit_hsps"),
        None,
    )
    hsps: List[BlastHSP] = []

    if hsp_container is not None:
        for hsp_node in hsp_container:
            if _local_tag(hsp_node) == "Hsp":
                hsps.append(_parse_hsp(hsp_node))

    definition = _child_text(hit_element, "Hit_def")
    scientific_names = []
    for tag in ("Hit_sciname", "Hit_scientific_name", "Hit_taxname"):
        value = _child_text(hit_element, tag)
        if value and value not in scientific_names:
            scientific_names.append(value)

    taxonomy_ids = []
    for tag in ("Hit_taxid", "Hit_taxid_list"):
        value = _child_text(hit_element, tag)
        if value:
            taxonomy_ids.extend(
                item.strip() for item in value.split(";") if item.strip()
            )

    # Legacy/compact NCBI XML may omit explicit taxonomy fields. Preserve the
    # real BLAST definition and add only a conservative, clearly documented
    # fallback label so the taxonomy overview remains useful.
    if not scientific_names:
        inferred_name = _infer_scientific_name_from_definition(definition)
        if inferred_name:
            scientific_names.append(inferred_name)

    return BlastHit(
        hit_id=_child_text(hit_element, "Hit_id"),
        accession=_child_text(hit_element, "Hit_accession"),
        definition=definition,
        length=_to_int(_child_text(hit_element, "Hit_len")),
        hsps=hsps,
        scientific_names=scientific_names,
        taxonomy_ids=taxonomy_ids,
    )

def _parse_parameters(root: ET.Element) -> Dict[str, Any]:
    parameters: Dict[str, Any] = {}
    for container in _descendants(root, "Parameters"):
        for child in container:
            tag = _local_tag(child)
            key = tag.replace("Parameters_", "")
            value = (child.text or "").strip()
            parameters[key] = value
        break
    return parameters


def _parse_database_statistics(root: ET.Element) -> Dict[str, Any]:
    stats: Dict[str, Any] = {}
    for container in _descendants(root, "Statistics"):
        for child in container:
            tag = _local_tag(child).replace("Statistics_", "")
            value = (child.text or "").strip()
            parsed: Any = value
            if re.fullmatch(r"-?\d+", value):
                parsed = int(value)
            elif re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)(?:[Ee][+-]?\d+)?", value):
                parsed = float(value)
            stats[tag] = parsed
        break
    return stats


def parse_blast_xml(
    xml_text: str,
    *,
    query: BlastQuery,
    rid: str,
    parameters: BlastParameters,
    submitted_at: str,
    completed_at: Optional[str] = None,
    database_type_override: Optional[str] = None,
) -> BlastRun:
    """Parse NCBI BLAST XML into a Synapse BlastRun."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise BlastAnalysisError(f"Invalid BLAST XML: {exc}") from exc

    program = (
        _child_text(root, "BlastOutput_program")
        or parameters.program
    )
    version = _child_text(root, "BlastOutput_version") or None
    database_reported = _child_text(root, "BlastOutput_db") or None

    hits: List[BlastHit] = []

    # Legacy XML has Iteration -> Iteration_hits -> Hit.
    # Search generically to remain tolerant of wrapper/container differences.
    for hit_element in _descendants(root, "Hit"):
        hits.append(_parse_hit(hit_element))

    metadata = BlastRunMetadata(
        rid=rid,
        submitted_at_utc=submitted_at,
        completed_at_utc=completed_at,
        program=program,
        database=parameters.database,
        database_type=(
            database_type_override
            or DATABASES.get(parameters.database, {}).get("type", "unknown")
        ),
        query_name=query.name,
        query_length=query.length,
        query_source=query.source,
        parameters=asdict(parameters),
        ncbi_status="READY",
        ncbi_version=version,
        ncbi_database_reported=database_reported,
    )

    notes = [
        "BLAST reports local sequence similarity and statistical significance.",
        "A high-scoring hit does not by itself prove biological identity or function.",
        "Interpretation depends on database composition, search parameters and query context.",
    ]

    if not hits:
        notes.append("No BLAST hits were returned for the submitted search.")

    return BlastRun(
        metadata=metadata,
        query=query,
        hits=hits,
        database_statistics=_parse_database_statistics(root),
        interpretation_notes=notes,
    )


# ============================================================================
# HIT ANALYSIS
# ============================================================================

def calculate_evidence_score(hit: BlastHit) -> HitEvidence:
    """
    Transparent heuristic evidence profile.

    This is NOT a BLAST score and must not be interpreted as experimental
    confidence. It summarizes returned alignment evidence for Synapse display.
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

    identity_component = min(100.0, max(0.0, hsp.identity_percent))
    coverage_component = min(100.0, max(0.0, hsp.query_coverage))

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

    if score >= 90:
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
    ]

    return HitEvidence(
        hit_id=hit.hit_id,
        accession=hit.accession,
        identity_percent=round(hsp.identity_percent, 3),
        query_coverage_percent=round(hsp.query_coverage, 3),
        evalue=hsp.evalue,
        bit_score=hsp.bit_score,
        alignment_length=hsp.alignment_length,
        evidence_score=score,
        evidence_band=band,
        notes=notes,
    )


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


def build_evidence_profile(run: BlastRun) -> List[HitEvidence]:
    return [calculate_evidence_score(hit) for hit in rank_hits(run.hits)]


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


# ============================================================================
# ORF / SYNAPSE INTEGRATION
# ============================================================================

def extract_orf_context(orf: Any) -> Optional[Tuple[str, int, str]]:
    """Adapt common Synapse ORF dict/object representations."""
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

    sequence = get_value(("sequence", "dna", "orf_sequence", "seq"))
    if not sequence:
        return None

    try:
        sequence = normalize_sequence(str(sequence), "DNA")
    except BlastAnalysisError:
        return None

    start = get_value(
        ("start", "start_position", "genomic_start", "position"),
        1,
    )
    try:
        start = int(start)
    except (TypeError, ValueError):
        start = 1

    name = str(
        get_value(("name", "id", "orf_id", "label"), "Selected ORF")
    )

    return sequence, max(1, start), name


def build_orf_query(
    orf: Any,
    *,
    sequence_name: str = "Selected Sequence",
) -> BlastQuery:
    extracted = extract_orf_context(orf)
    if extracted is None:
        raise BlastAnalysisError("The selected ORF could not be extracted.")

    sequence, start, name = extracted
    return build_query(
        sequence,
        name=name,
        sequence_type="DNA",
        source="Selected ORF",
        region_start=start,
        region_end=start + len(sequence) - 1,
        feature_type="ORF",
        feature_name=name,
    )


def run_blast_analysis(
    sequence: str,
    *,
    sequence_name: str = "Selected Sequence",
    sequence_type: Optional[str] = None,
    source: str = "Selected Sequence",
    selected_orf: Optional[Any] = None,
    program: str = "blastn",
    database: str = "nt",
    evalue: float = 10.0,
    max_target_seqs: int = DEFAULT_MAX_HITS,
    email: str = DEFAULT_EMAIL,
    tool: str = DEFAULT_TOOL,
    engine: str = "remote",
    local_database_path: Optional[str | Path] = None,
    local_database_root: Optional[str | Path] = None,
    client: Optional[NCBIBlastClient] = None,
) -> BlastRun:
    """
    Programmatic Synapse entry point.

    If selected_orf is supplied, the ORF becomes the BLAST query while global
    coordinates are retained in the query metadata.

    engine:
        "remote" -> NCBI BLAST URL API
        "local"  -> installed NCBI BLAST+ against a local database
    """
    if selected_orf is not None:
        query = build_orf_query(
            selected_orf,
            sequence_name=sequence_name,
        )
    else:
        query = build_query(
            sequence,
            name=sequence_name,
            sequence_type=sequence_type,
            source=source,
        )

    parameters = BlastParameters(
        program=program,
        database=database,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        email=email,
        tool=tool,
    )

    engine_normalized = engine.lower().strip()
    if engine_normalized == "remote":
        parameters = parameters.normalized()
        api_client = client or NCBIBlastClient(
            tool=parameters.tool,
            email=parameters.email,
        )
        return api_client.run(query, parameters)

    if engine_normalized != "local":
        raise BlastAnalysisError("engine must be 'local' or 'remote'.")

    manager = LocalBlastDatabaseManager(root=local_database_root)
    if local_database_path:
        db_path = Path(local_database_path).expanduser().resolve()
        if not db_path.exists():
            raise BlastAnalysisError(f"Local BLAST database path not found: {db_path}")
        database_type = (
            LOCAL_DATABASE_CATALOG.get(database, {}).get("type")
            or "nucleotide"
        )
        database_record = LocalBlastDatabase(
            name=database,
            database_type=database_type,
            path=str(db_path),
            installed=True,
            structurally_valid=True,
            verification_status="USER_SUPPLIED",
            verification_message="Database path supplied programmatically.",
        )
    else:
        database_record = manager.scan_database(database)

    if not database_record.exists:
        raise BlastAnalysisError(
            f"Local BLAST database '{database}' is not available under "
            f"{manager.root}."
        )

    if database_record.database_type not in {"nucleotide", "protein"}:
        raise BlastAnalysisError("Local database type is invalid.")

    expected_db_type = PROGRAMS.get(program.lower().strip(), {}).get("database_type")
    if expected_db_type != database_record.database_type:
        raise BlastAnalysisError(
            f"{program} requires a {expected_db_type} database, "
            f"but {database_record.name} is {database_record.database_type}."
        )

    # Local validation is intentionally independent from the remote catalog.
    if parameters.evalue <= 0 or parameters.max_target_seqs < 1:
        raise BlastAnalysisError("Invalid BLAST parameters.")

    engine_client = LocalBlastEngine()
    return engine_client.run(query, parameters, database_record)


# ============================================================================
# SERIALIZATION / EXPORT
# ============================================================================

def hsp_to_dict(hsp: BlastHSP) -> Dict[str, Any]:
    data = asdict(hsp)
    data["identity_percent"] = hsp.identity_percent
    data["positive_percent"] = hsp.positive_percent
    data["gap_percent"] = hsp.gap_percent
    data["query_coverage"] = hsp.query_coverage
    return data


def hit_to_dict(hit: BlastHit) -> Dict[str, Any]:
    data = asdict(hit)
    data["best_bit_score"] = hit.best_bit_score
    data["best_evalue"] = hit.best_evalue
    data["best_identity_percent"] = hit.best_identity_percent
    data["best_coverage_percent"] = hit.best_coverage_percent
    data["total_alignment_length"] = hit.total_alignment_length
    data["hsps"] = [hsp_to_dict(hsp) for hsp in hit.hsps]
    return data


def run_to_dict(run: BlastRun) -> Dict[str, Any]:
    return {
        "metadata": asdict(run.metadata),
        "query": asdict(run.query),
        "hits": [hit_to_dict(hit) for hit in rank_hits(run.hits)],
        "database_statistics": run.database_statistics,
        "interpretation_notes": run.interpretation_notes,
        "evidence_profile": [
            asdict(profile) for profile in build_evidence_profile(run)
        ],
    }


def run_to_json(run: BlastRun, *, indent: int = 2) -> str:
    return json.dumps(
        run_to_dict(run),
        indent=indent,
        ensure_ascii=False,
    )


def export_json(run: BlastRun, path: str | Path) -> str:
    path = _resolve_synapse_user_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(run_to_json(run), encoding="utf-8")
    return str(path)


def export_csv(run: BlastRun, path: str | Path) -> str:
    path = _resolve_synapse_user_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "rank",
        "accession",
        "hit_id",
        "definition",
        "identity_percent",
        "query_coverage_percent",
        "evalue",
        "bit_score",
        "alignment_length",
        "evidence_score",
        "evidence_band",
    ]

    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        profiles = {
            profile.accession: profile
            for profile in build_evidence_profile(run)
        }

        for rank, hit in enumerate(rank_hits(run.hits), start=1):
            profile = profiles.get(hit.accession)
            writer.writerow(
                {
                    "rank": rank,
                    "accession": hit.accession,
                    "hit_id": hit.hit_id,
                    "definition": hit.definition,
                    "identity_percent": (
                        profile.identity_percent if profile else hit.best_identity_percent
                    ),
                    "query_coverage_percent": (
                        profile.query_coverage_percent if profile else hit.best_coverage_percent
                    ),
                    "evalue": hit.best_evalue,
                    "bit_score": hit.best_bit_score,
                    "alignment_length": (
                        profile.alignment_length
                        if profile
                        else hit.best_hsp.alignment_length if hit.best_hsp else 0
                    ),
                    "evidence_score": profile.evidence_score if profile else "",
                    "evidence_band": profile.evidence_band if profile else "",
                }
            )

    return str(path)


def build_text_report(run: BlastRun) -> str:
    lines = [
        "=" * 90,
        "SYNAPSE BLAST REPORT",
        "=" * 90,
        f"Query:            {run.query.name}",
        f"Query type:       {run.query.sequence_type}",
        f"Query length:     {run.query.length}",
        f"Query source:     {run.query.source}",
        f"Feature type:     {run.query.feature_type}",
        f"Program:          {run.metadata.program}",
        f"Database:          {run.metadata.database}",
        f"RID:              {run.metadata.rid}",
        f"Hits:             {run.hit_count}",
        f"Significant hits: {run.significant_hit_count}",
        "",
        "TOP HITS",
        "-" * 90,
    ]

    for index, hit in enumerate(rank_hits(run.hits)[:20], start=1):
        lines.extend(
            [
                f"[{index}] {hit.accession or hit.hit_id}",
                f"    Definition: {hit.definition}",
                f"    Identity:   {hit.best_identity_percent:.2f} %",
                f"    Coverage:   {hit.best_coverage_percent:.2f} %",
                f"    E-value:    {hit.best_evalue}",
                f"    Bit score:  {hit.best_bit_score:.2f}",
                "",
            ]
        )

    lines.extend(
        [
            "SCIENTIFIC INTERPRETATION NOTES",
            "-" * 90,
            *[f"- {note}" for note in run.interpretation_notes],
            "",
            "RUN METADATA",
            "-" * 90,
            json.dumps(asdict(run.metadata), indent=2, ensure_ascii=False),
        ]
    )

    return "\r\n".join(lines)


def export_text(run: BlastRun, path: str | Path) -> str:
    path = _resolve_synapse_user_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_text_report(run), encoding="utf-8")
    return str(path)


# ============================================================================
# RUN COMPARISON
# ============================================================================

def compare_blast_runs(run_a: BlastRun, run_b: BlastRun) -> BlastComparison:
    map_a = {
        hit.accession or hit.hit_id: hit
        for hit in run_a.hits
    }
    map_b = {
        hit.accession or hit.hit_id: hit
        for hit in run_b.hits
    }

    keys_a = set(map_a)
    keys_b = set(map_b)
    common = sorted(keys_a & keys_b)

    identity_differences: Dict[str, float] = {}
    coverage_differences: Dict[str, float] = {}

    for accession in common:
        identity_differences[accession] = round(
            map_a[accession].best_identity_percent
            - map_b[accession].best_identity_percent,
            4,
        )
        coverage_differences[accession] = round(
            map_a[accession].best_coverage_percent
            - map_b[accession].best_coverage_percent,
            4,
        )

    return BlastComparison(
        query_a=run_a.query.name,
        query_b=run_b.query.name,
        common_accessions=common,
        unique_to_a=sorted(keys_a - keys_b),
        unique_to_b=sorted(keys_b - keys_a),
        identity_differences=identity_differences,
        coverage_differences=coverage_differences,
    )


# ============================================================================
# DISPLAY
# ============================================================================

def _line(char: str = "-", width: int = 90) -> str:
    return char * width


def print_hit_summary(run: BlastRun, *, limit: int = 15) -> None:
    print("\r\n" + _line("="))
    print("BLAST RESULT SUMMARY")
    print(_line("="))
    print(f"Query:              {run.query.name}")
    print(f"Type:               {run.query.sequence_type}")
    print(f"Length:             {run.query.length}")
    print(f"Program:            {run.metadata.program}")
    print(f"Database:            {run.metadata.database}")
    print(f"RID:                {run.metadata.rid}")
    print(f"Hits:               {run.hit_count}")
    print(f"Significant hits:   {run.significant_hit_count}")
    print(_line())

    print(
        f"{'Rank':<6}{'Accession':<18}{'Identity':>12}"
        f"{'Coverage':>12}{'E-value':>16}{'Bit score':>14}"
    )
    print(_line())

    for rank, hit in enumerate(rank_hits(run.hits)[:limit], start=1):
        evalue = f"{hit.best_evalue:.3g}" if hit.best_evalue is not None else "n/a"
        print(
            f"{rank:<6}"
            f"{(hit.accession or hit.hit_id)[:17]:<18}"
            f"{hit.best_identity_percent:>11.2f}%"
            f"{hit.best_coverage_percent:>11.2f}%"
            f"{evalue:>16}"
            f"{hit.best_bit_score:>14.2f}"
        )


def print_hit_detail(hit: BlastHit) -> None:
    print("\r\n" + _line("="))
    print("BLAST HIT DETAIL")
    print(_line("="))
    print(f"Accession:      {hit.accession}")
    print(f"Hit ID:         {hit.hit_id}")
    print(f"Definition:     {hit.definition}")
    print(f"Subject length: {hit.length}")
    print(f"HSPs:           {len(hit.hsps)}")
    print()

    for hsp in hit.hsps:
        print(
            f"HSP {hsp.number}: identity={hsp.identity_percent:.2f}% | "
            f"coverage={hsp.query_coverage:.2f}% | "
            f"E-value={hsp.evalue:.3g} | "
            f"bit score={hsp.bit_score:.2f}"
        )


# ============================================================================
# INTERACTIVE MENU
# ============================================================================

def _print_local_database_status(manager: LocalBlastDatabaseManager) -> None:
    print("\r\nLOCAL BLAST DATABASE STATUS")
    print(_line("="))
    print(f"Root: {manager.root}")
    print(_line())
    print(
        f"{'Database':<24}{'Type':<12}{'Status':<42}"
        f"{'Volumes':>12}{'Size':>12}"
    )
    print(_line())

    for database in manager.list_databases():
        if database.exists:
            status = database.verification_status
            size = (
                f"{database.size_gb:.2f} GB"
                if database.size_gb >= 1
                else f"{database.size_mb:.1f} MB"
            )
        elif database.installed:
            if database.missing_volume_files:
                first_volume, first_missing = next(iter(database.missing_volume_files.items()))
                status = f"INCOMPLETE ({first_volume}: {','.join(first_missing)})"
            elif database.missing_volumes:
                status = f"INCOMPLETE (missing {len(database.missing_volumes)})"
            else:
                status = "INCOMPLETE"
            size = (
                f"{database.size_gb:.2f} GB"
                if database.size_gb >= 1
                else f"{database.size_mb:.1f} MB"
            )
        else:
            status = "NOT INSTALLED"
            size = "-"

        if database.expected_volume_count is not None:
            volume_text = (
                f"{database.volume_count}/{database.expected_volume_count}"
            )
        elif database.volume_count:
            volume_text = str(database.volume_count)
        else:
            volume_text = "1" if database.installed else "-"

        print(
            f"{database.name:<24}"
            f"{database.database_type:<12}"
            f"{status:<24}"
            f"{volume_text:>12}"
            f"{size:>12}"
        )


def _ensure_local_database_verified(
    manager: LocalBlastDatabaseManager,
    database: LocalBlastDatabase,
) -> LocalBlastDatabase:
    """Ensure a selected local database passes a real BLASTDB verification.

    Structural discovery alone is not enough for a multi-volume NCBI database.
    If the database has not yet been verified, run ``blastdbcmd -info`` through
    the same runtime-view mechanism used by the database manager. This prevents
    a search from starting against a database that Synapse can see on disk but
    BLAST+ cannot actually open.
    """
    if database.verification_status == "VERIFIED":
        return database

    verified = manager.verify_database(database.name, run_blastdbcmd=True)
    if verified.verification_status != "VERIFIED":
        message = verified.verification_message or "No verification diagnostic was returned."
        raise BlastAnalysisError(
            f"Local BLAST database '{database.name}' is not verified. {message}"
        )
    return verified


def _select_local_database(
    manager: LocalBlastDatabaseManager,
    database_type: str,
) -> Optional[LocalBlastDatabase]:
    available = [
        database
        for database in manager.list_databases()
        if database.database_type == database_type and database.exists
    ]

    if not available:
        print(
            "\r\nNo compatible local BLAST database is installed.\r\n"
            "Use [3] BLAST Database Manager to download or create one."
        )
        return None

    # Prefer a fully verified database and, among equal states, prefer the
    # canonical NCBI RefSeq database over small test/custom databases.
    status_rank = {
        "VERIFIED": 0,
        "STRUCTURALLY_VALID": 1,
        "USER_SUPPLIED": 2,
    }
    available.sort(
        key=lambda database: (
            status_rank.get(database.verification_status, 3),
            0 if database.name.lower() == "refseq_rna" else 1,
            -database.size_bytes,
            database.name.lower(),
        )
    )

    print("\r\nLOCAL DATABASE")
    for index, database in enumerate(available, start=1):
        if database.verification_status == "VERIFIED":
            status_label = "VERIFIED"
        else:
            status_label = database.verification_status
        sequence_info = (
            f" | {database.sequence_count:,} seq"
            if database.sequence_count is not None
            else ""
        )
        print(
            f"[{index}] {database.name} | "
            f"{database.size_gb:.2f} GB | {status_label}{sequence_info}"
        )
    print("[0] Cancel")

    choice = input("Database: ").strip()
    if choice == "0":
        return None

    try:
        return available[int(choice) - 1]
    except (ValueError, IndexError):
        print("Invalid database selection.")
        return None


def blast_database_manager_menu(
    manager: Optional[LocalBlastDatabaseManager] = None,
) -> None:
    """Interactive manager for local BLAST+ databases."""
    manager = manager or LocalBlastDatabaseManager()

    while True:
        print("\r\n" + _line("="))
        print("BLAST DATABASE MANAGER")
        print(_line("="))
        print(f"Database root: {manager.root}")
        print()
        print("DATABASE")
        print("[1] Database Status")
        print("[2] Download / Update NCBI Database")
        print("[3] Verify Database")
        print("[4] Create Custom Database from FASTA")
        print("[5] Remove Database")
        print("[6] Database Information")
        print()
        print("TOOLS")
        print("[7] BLAST+ Tool Status")
        print("[8] Change Database Root")
        print("[0] Back")
        print(_line())

        choice = input("Selection: ").strip()

        if choice == "0":
            return

        if choice == "1":
            _print_local_database_status(manager)
            continue

        if choice == "2":
            print("\r\nNCBI DATABASES AVAILABLE THROUGH update_blastdb.pl")
            print(_line())
            for index, (name, info) in enumerate(
                LOCAL_DATABASE_CATALOG.items(), start=1
            ):
                recommendation = "RECOMMENDED START" if info["recommended"] else "ADVANCED"
                print(
                    f"[{index}] {name:<18} "
                    f"{info['type']:<12} "
                    f"{recommendation:<16} "
                    f"{info['description']}"
                )
            print("[0] Cancel")

            raw = input("Database: ").strip()
            if raw == "0":
                continue

            names = list(LOCAL_DATABASE_CATALOG)
            try:
                name = names[int(raw) - 1]
            except (ValueError, IndexError):
                print("Invalid database selection.")
                continue

            print()
            print(f"Selected database: {name}")
            print(
                "IMPORTANT: The database will be downloaded from NCBI into the "
                f"following directory:\r\n{manager.root}"
            )
            print(
                "Large databases can require substantial disk space and network "
                "bandwidth. Synapse does not start this automatically."
            )
            confirm = input("Type DOWNLOAD to continue: ").strip()
            if confirm != "DOWNLOAD":
                print("Download cancelled.")
                continue

            try:
                database = manager.update_ncbi_database(name)
                print(f"\r\nUpdate completed: {database.name}")
                print(f"Status: {database.verification_status}")
                print(f"Files:  {database.file_count}")
                print(f"Size:   {database.size_gb:.3f} GB")
            except BlastAnalysisError as exc:
                print(f"ERROR: {exc}")
            continue

        if choice == "3":
            _print_local_database_status(manager)
            name = input("\r\nDatabase name to verify (0 cancel): ").strip()
            if name == "0":
                continue
            try:
                database = manager.verify_database(name, run_blastdbcmd=True)
                print("\r\nVERIFICATION RESULT")
                print(_line())
                print(f"Database: {database.name}")
                print(f"Status:   {database.verification_status}")
                print(f"Message:  {database.verification_message}")
            except BlastAnalysisError as exc:
                print(f"ERROR: {exc}")
            continue

        if choice == "4":
            fasta = input("FASTA file path: ").strip()
            if not fasta:
                print("No FASTA path supplied.")
                continue

            name = input("New database name: ").strip()
            if not name:
                print("No database name supplied.")
                continue

            print("\r\nDATABASE TYPE")
            print("[1] Nucleotide")
            print("[2] Protein")
            print("[0] Cancel")
            db_choice = input("Type: ").strip()
            if db_choice == "0":
                continue
            database_type = {"1": "nucleotide", "2": "protein"}.get(db_choice)
            if database_type is None:
                print("Invalid database type.")
                continue

            try:
                database = manager.create_custom_database(
                    fasta,
                    name,
                    database_type,
                )
                print("\r\nCUSTOM DATABASE CREATED")
                print(_line())
                print(f"Name:   {database.name}")
                print(f"Type:   {database.database_type}")
                print(f"Path:   {database.path}")
                print(f"Files:  {database.file_count}")
                print(f"Size:   {database.size_mb:.2f} MB")
            except BlastAnalysisError as exc:
                print(f"ERROR: {exc}")
            continue

        if choice == "5":
            _print_local_database_status(manager)
            name = input("\r\nDatabase name to remove (0 cancel): ").strip()
            if name == "0":
                continue

            database = manager.scan_database(name)
            if not database.installed:
                print("Database is not installed.")
                continue

            confirm = input(
                f"Type REMOVE {database.name} to permanently remove its local files: "
            ).strip()
            if confirm != f"REMOVE {database.name}":
                print("Removal cancelled.")
                continue

            try:
                removed = manager.remove_database(name)
                print(f"Removed {len(removed)} files.")
            except BlastAnalysisError as exc:
                print(f"ERROR: {exc}")
            continue

        if choice == "6":
            _print_local_database_status(manager)
            name = input("\r\nDatabase name (0 cancel): ").strip()
            if name == "0":
                continue
            try:
                info = manager.database_info(name)
                print("\r\nDATABASE INFORMATION")
                print(_line())
                print(json.dumps(info, indent=2, ensure_ascii=False))
            except BlastAnalysisError as exc:
                print(f"ERROR: {exc}")
            continue

        if choice == "7":
            print("\r\nBLAST+ TOOL STATUS")
            print(_line())
            tools = manager.executable_manager.available()
            for name, path in tools.items():
                print(f"{name:<14}: {path or 'NOT FOUND'}")
            print(
                f"perl           : "
                f"{manager.executable_manager.resolve_perl() or 'NOT FOUND'}"
            )
            print(
                f"{BLAST_UPDATE_SCRIPT:<14}: "
                f"{manager.executable_manager.resolve_update_script() or 'NOT FOUND'}"
            )
            continue

        if choice == "8":
            raw_root = input(
                "New database root (blank = default project data/blast_db, 0 cancel): "
            ).strip()
            if raw_root == "0":
                continue
            try:
                manager = LocalBlastDatabaseManager(root=raw_root or DEFAULT_LOCAL_DB_ROOT)
                manager.ensure_root()
                print(f"Database root changed to: {manager.root}")
                print(
                    f"To make this persistent across sessions, set the "
                    f"{LOCAL_DB_ROOT_ENV} environment variable."
                )
            except OSError as exc:
                print(f"ERROR: {exc}")
            continue

        print("Invalid selection.")


def blast_menu(
    sequence: Optional[str] = None,
    sequence_name: str = "Selected Sequence",
    selected_orf: Optional[Any] = None,
) -> None:
    """
    Interactive Synapse BLAST Search menu.

    Synapse supports two explicit search engines:
        1. Local BLAST+ against a managed local database.
        2. NCBI remote BLAST through the NCBI URL API.

    If an ORF is selected, it is used as the BLAST query context and its
    global coordinates are retained. A user can also load an external FASTA
    file as an explicit query source; no implicit first-record-only behavior
    is used for multi-record FASTA files.
    """
    base_sequence = sequence
    current_name = sequence_name
    current_orf = selected_orf
    active_query: Optional[BlastQuery] = None
    last_run: Optional[BlastRun] = None
    comparison_run: Optional[BlastRun] = None
    manager = LocalBlastDatabaseManager()

    while True:
        print("\r\n" + _line("="))
        print("BLAST SEARCH")
        print(_line("="))

        if active_query is not None:
            print(f"Query context:    {active_query.name}")
            print(f"Query source:     {active_query.source}")
            print(f"Query length:     {active_query.length} {active_query.sequence_type.lower()}")
        elif current_orf is not None:
            try:
                orf_query = build_orf_query(current_orf, sequence_name=current_name)
                print(f"Query context:    Selected ORF ({orf_query.name})")
                print(f"Query length:     {orf_query.length} bp")
            except BlastAnalysisError:
                print("Query context:    Selected ORF (unavailable)")
        elif base_sequence:
            print(f"Query context:    {current_name}")
            print("Query source:     Selected Synapse Sequence")
            print(f"Query length:     {len(''.join(base_sequence.split()))} nt")
        else:
            print("Query context:    No sequence loaded")

        if last_run:
            print(f"Last engine:      {last_run.metadata.engine}")
            print(f"Last database:    {last_run.metadata.database}")

        print()
        print("SEARCH")
        print("[1] Run BLAST Search")
        print("[2] Search Configuration")
        print("[3] BLAST Database Manager")
        print()
        print("RESULT ANALYSIS")
        print("[4] Result Summary")
        print("[5] Hit Detail & Alignment")
        print("[6] Evidence Profile")
        print("[7] Taxonomic Overview")
        print("[8] Compare BLAST Runs")
        print()
        print("OUTPUT")
        print("[9] Export Results")
        print("[10] Full BLAST Report")
        print()
        print("SYSTEM")
        print("[11] BLAST Self-Test")
        print("[12] Supported Programs & Databases")
        print("[13] Database Status")
        print("[0] Back")
        print(_line())

        choice = input("Selection: ").strip()

        if choice == "0":
            return

        if choice == "1":
            print("\r\nQUERY SOURCE")
            print("[1] Selected Synapse Sequence")
            print("[2] Load FASTA File")
            print("[0] Cancel")
            query_source_choice = input("Source: ").strip()

            if query_source_choice == "0":
                continue

            search_query: Optional[BlastQuery] = None

            if query_source_choice == "1":
                active_query = None
                if not base_sequence and current_orf is None:
                    print("ERROR: No DNA/protein sequence or ORF is available.")
                    continue

                try:
                    if current_orf is not None:
                        search_query = build_orf_query(
                            current_orf,
                            sequence_name=current_name,
                        )
                    else:
                        selected_type = infer_sequence_type(base_sequence or "")
                        search_query = build_query(
                            base_sequence or "",
                            name=current_name,
                            sequence_type=selected_type,
                            source="Selected Synapse Sequence",
                        )
                except BlastAnalysisError as exc:
                    print(f"ERROR: {exc}")
                    continue

            elif query_source_choice == "2":
                raw_path = input("FASTA file path: ").strip().strip('\"')
                if not raw_path:
                    print("No FASTA path supplied.")
                    continue

                try:
                    fasta_queries = load_fasta_queries(raw_path)
                except BlastAnalysisError as exc:
                    print(f"ERROR: {exc}")
                    continue

                if len(fasta_queries) == 1:
                    search_query = fasta_queries[0]
                else:
                    print("\r\nFASTA RECORDS")
                    for index, fasta_query in enumerate(fasta_queries, start=1):
                        print(
                            f"[{index}] {fasta_query.name} | "
                            f"{fasta_query.length} {fasta_query.sequence_type.lower()}"
                        )
                    print("[0] Cancel")
                    record_choice = input("Record: ").strip()
                    if record_choice == "0":
                        continue
                    try:
                        search_query = fasta_queries[int(record_choice) - 1]
                    except (ValueError, IndexError):
                        print("Invalid FASTA record.")
                        continue

                active_query = search_query
                print(
                    f"Loaded FASTA query: {search_query.name} "
                    f"({search_query.length} {search_query.sequence_type.lower()})"
                )

            else:
                print("Invalid query source.")
                continue

            if search_query is None:
                print("ERROR: No BLAST query was selected.")
                continue

            available_query_type = search_query.sequence_type

            compatible_programs = [
                name
                for name, info in PROGRAMS.items()
                if info["query_type"] == available_query_type
            ]

            print("\r\nSEARCH ENGINE")
            print("[1] Local BLAST+")
            print("[2] NCBI BLAST")
            print("[0] Cancel")
            engine_choice = input("Engine: ").strip()
            if engine_choice == "0":
                continue

            engine = {"1": "local", "2": "remote"}.get(engine_choice)
            if engine is None:
                print("Invalid engine.")
                continue

            print("\r\nBLAST PROGRAM")
            for index, name in enumerate(compatible_programs, start=1):
                print(
                    f"[{index}] {name.upper()} - "
                    f"{PROGRAMS[name]['description']}"
                )
            print("[0] Cancel")

            program_choice = input("Program: ").strip()
            if program_choice == "0":
                continue

            try:
                program_index = int(program_choice)
                program = compatible_programs[program_index - 1]
            except (ValueError, IndexError):
                print("Invalid program.")
                continue

            database: Optional[str] = None
            local_database: Optional[LocalBlastDatabase] = None

            if engine == "local":
                database_type = PROGRAMS[program]["database_type"]
                local_database = _select_local_database(manager, database_type)
                if local_database is None:
                    continue
                try:
                    local_database = _ensure_local_database_verified(
                        manager, local_database
                    )
                except BlastAnalysisError as exc:
                    print(f"ERROR: {exc}")
                    continue
                database = local_database.name
                print(
                    f"Verified local database: {local_database.name} | "
                    f"{local_database.sequence_count:,} sequences"
                    if local_database.sequence_count is not None
                    else f"Verified local database: {local_database.name}"
                )
            else:
                compatible_databases = [
                    name
                    for name, info in DATABASES.items()
                    if info["type"] == PROGRAMS[program]["database_type"]
                ]

                print("\r\nNCBI REMOTE DATABASE")
                for index, name in enumerate(compatible_databases, start=1):
                    print(
                        f"[{index}] {name} - "
                        f"{DATABASES[name]['description']}"
                    )
                print("[0] Cancel")

                db_choice = input("Database: ").strip()
                if db_choice == "0":
                    continue

                try:
                    db_index = int(db_choice)
                    database = compatible_databases[db_index - 1]
                except (ValueError, IndexError):
                    print("Invalid database.")
                    continue

            raw_evalue = input("E-value [10]: ").strip()
            raw_max = input(f"Maximum target sequences [{DEFAULT_MAX_HITS}]: ").strip()

            raw_threads = ""
            raw_word_size = ""
            raw_filter = ""
            if engine == "local":
                raw_threads = input("CPU threads [4]: ").strip()
                raw_word_size = input("Word size [default]: ").strip()
                raw_filter = input(
                    "Filtering [default/none]: "
                ).strip().lower()

            try:
                evalue = float(raw_evalue) if raw_evalue else 10.0
                max_hits = int(raw_max) if raw_max else DEFAULT_MAX_HITS
                threads = int(raw_threads) if raw_threads else 4
                word_size = int(raw_word_size) if raw_word_size else None
                filter_query = None
                if raw_filter in {"none", "no", "off"}:
                    filter_query = "no"
                elif raw_filter in {"default", "", "yes", "on"}:
                    filter_query = None
                else:
                    raise BlastAnalysisError(
                        "Filtering must be 'default' or 'none'."
                    )
                if threads < 1:
                    raise BlastAnalysisError("CPU threads must be at least 1.")
                if word_size is not None and word_size < 1:
                    raise BlastAnalysisError("Word size must be at least 1.")

                query = search_query
                if query is None:
                    raise BlastAnalysisError("No BLAST query was selected.")

                parameters = BlastParameters(
                    program=program,
                    database=database or "",
                    evalue=evalue,
                    max_target_seqs=max_hits,
                    word_size=word_size,
                    filter_query=filter_query,
                    num_threads=threads if engine == "local" else None,
                )

                if engine == "local":
                    assert local_database is not None
                    print("\r\nRunning local BLAST+ search...")
                    new_run = LocalBlastEngine().run(
                        query,
                        parameters,
                        local_database,
                    )
                else:
                    normalized = parameters.normalized()
                    print("\r\nSubmitting BLAST search to NCBI...")
                    print("The search may take several seconds or longer.")
                    client = NCBIBlastClient(
                        tool=normalized.tool,
                        email=normalized.email,
                    )
                    new_run = client.run(query, normalized)

                comparison_run = last_run
                last_run = new_run
                print("\r\nBLAST SEARCH COMPLETED")
                print_hit_summary(last_run)
                if not last_run.hits:
                    print(
                        "\r\nNo hits were returned. This is a valid BLAST result; "
                        "it does not indicate that the query is biologically irrelevant."
                    )

            except (BlastAnalysisError, AssertionError) as exc:
                print(f"ERROR: {exc}")

        elif choice == "2":
            print("\r\nSEARCH CONFIGURATION")
            print(_line())
            print("Supported programs:")
            for name, info in PROGRAMS.items():
                print(
                    f"- {name.upper():<8} "
                    f"query={info['query_type']:<8} "
                    f"database={info['database_type']:<11}"
                )
            print()
            print(f"Default E-value:           {10.0}")
            print(f"Default max target seqs:   {DEFAULT_MAX_HITS}")
            print("Default local CPU threads: 4")
            print("Local filtering:           BLAST+ default")
            print(f"Local database root:        {manager.root}")
            print(
                "NCBI remote searches use the NCBI URL API. "
                "Local searches use installed BLAST+ executables."
            )

        elif choice == "3":
            blast_database_manager_menu(manager)
            manager = LocalBlastDatabaseManager(root=manager.root)

        elif choice == "4":
            if not last_run:
                print("Run [1] BLAST Search first.")
                continue
            print_hit_summary(last_run)

        elif choice == "5":
            if not last_run or not last_run.hits:
                print("No BLAST hits available.")
                continue

            ranked = rank_hits(last_run.hits)
            for index, hit in enumerate(ranked[:20], start=1):
                print(
                    f"[{index}] {hit.accession or hit.hit_id} | "
                    f"{hit.best_identity_percent:.2f}% identity | "
                    f"{hit.best_coverage_percent:.2f}% coverage"
                )

            raw = input("Hit number (0 cancel): ").strip()
            if raw == "0":
                continue

            try:
                index = int(raw)
                hit = ranked[index - 1]
                print_hit_detail(hit)
                print()
                print(alignment_text(hit))
            except (ValueError, IndexError, BlastAnalysisError) as exc:
                print(f"ERROR: {exc}")

        elif choice == "6":
            if not last_run:
                print("Run [1] BLAST Search first.")
                continue

            print("\r\nEVIDENCE PROFILE")
            print(_line())
            for index, profile in enumerate(build_evidence_profile(last_run), start=1):
                print(
                    f"[{index}] {profile.accession or profile.hit_id}: "
                    f"{profile.evidence_score:.2f}/100 | "
                    f"{profile.evidence_band}"
                )
                for note in profile.notes:
                    print(f"    - {note}")

        elif choice == "7":
            if not last_run:
                print("Run [1] BLAST Search first.")
                continue

            groups = group_hits_by_taxonomy(last_run.hits)
            print("\r\nTAXONOMIC OVERVIEW")
            print(_line())

            identified = sum(1 for hit in last_run.hits if hit.scientific_names)
            unknown = len(last_run.hits) - identified
            inferred = sum(
                1
                for hit in last_run.hits
                if hit.scientific_names
                and not hit.taxonomy_ids
                and _taxonomy_source(hit) == "inferred from BLAST hit definition"
            )

            print(f"Hits analyzed:       {len(last_run.hits)}")
            print(f"Taxa identified:     {identified}")
            print(f"Unknown:             {unknown}")
            if inferred:
                print(
                    f"Definition-inferred: {inferred} "
                    "(not an authoritative taxonomy assignment)"
                )
            print()

            if not groups:
                print("No taxonomy labels were available from this result set.")
            else:
                print("Taxon / organism                         Hits")
                print("-----------------------------------------------")
                for name, count in groups.items():
                    print(f"{name:<40} {count:>4}")

            if inferred:
                print()
                print(
                    "NOTE: Some organism names were inferred conservatively "
                    "from BLAST hit definitions because explicit taxonomy "
                    "fields were not present in the returned result."
                )

        elif choice == "8":
            if not last_run:
                print("Run [1] BLAST Search first.")
                continue

            if comparison_run is None:
                print(
                    "A second BLAST run is required. "
                    "Run another search first."
                )
                continue

            comparison = compare_blast_runs(comparison_run, last_run)
            print("\r\nBLAST RUN COMPARISON")
            print(_line())
            print(f"Query A:       {comparison.query_a}")
            print(f"Query B:       {comparison.query_b}")
            print(f"Common hits:   {len(comparison.common_accessions)}")
            print(f"Unique to A:   {len(comparison.unique_to_a)}")
            print(f"Unique to B:   {len(comparison.unique_to_b)}")

        elif choice == "9":
            if not last_run:
                print("Run [1] BLAST Search first.")
                continue

            path = input("Export path (.json/.csv/.txt): ").strip()
            if not path:
                print("No path supplied.")
                continue

            try:
                lower = path.lower()
                if lower.endswith(".json"):
                    export_json(last_run, path)
                elif lower.endswith(".csv"):
                    export_csv(last_run, path)
                elif lower.endswith(".txt"):
                    export_text(last_run, path)
                else:
                    print("Use .json, .csv or .txt.")
                    continue
                print(f"Exported: {path}")
            except OSError as exc:
                print(f"Export failed: {exc}")

        elif choice == "10":
            if not last_run:
                print("Run [1] BLAST Search first.")
            else:
                print(build_text_report(last_run))

        elif choice == "11":
            run_self_test(verbose=True)

        elif choice == "12":
            print("\r\nBLAST PROGRAMS")
            print(_line())
            for name, info in PROGRAMS.items():
                print(
                    f"{name.upper():<9} "
                    f"query={info['query_type']:<8} "
                    f"database={info['database_type']:<11} "
                    f"{info['description']}"
                )

            print("\r\nREMOTE DATABASES")
            print(_line())
            for name, info in DATABASES.items():
                print(f"{name:<18} {info['type']:<12} {info['description']}")

            print("\r\nLOCAL NCBI DATABASE CATALOG")
            print(_line())
            for name, info in LOCAL_DATABASE_CATALOG.items():
                print(
                    f"{name:<18} {info['type']:<12} "
                    f"{info['description']}"
                )

        elif choice == "13":
            _print_local_database_status(manager)

        else:
            print("Invalid selection.")


# ============================================================================
# SELF-TEST SUITE
# ============================================================================

SAMPLE_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.17.0+</BlastOutput_version>
  <BlastOutput_db>nt</BlastOutput_db>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>ref|TEST001|</Hit_id>
          <Hit_def>Example reference sequence</Hit_def>
          <Hit_accession>TEST001</Hit_accession>
          <Hit_len>120</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>100</Hsp_bit-score>
              <Hsp_score>110</Hsp_score>
              <Hsp_evalue>1e-25</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>20</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>29</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>19</Hsp_identity>
              <Hsp_positive>19</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>20</Hsp_align-len>
              <Hsp_qseq>ATGCGTACGTACGTACGTAC</Hsp_qseq>
              <Hsp_hseq>ATGCGTACGTACGTACGTAT</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||| </Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_statistics>
    <Statistics>
      <Statistics_db-num>1000</Statistics_db-num>
      <Statistics_db-len>1000000</Statistics_db-len>
    </Statistics>
  </BlastOutput_statistics>
</BlastOutput>
"""


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _test_sequence_validation() -> None:
    _assert(normalize_sequence(" a c g t n ") == "ACGTN", "DNA normalization failed")
    try:
        normalize_sequence("ACGTZ", "DNA")
    except BlastAnalysisError:
        pass
    else:
        raise AssertionError("Invalid DNA was accepted")


def _test_type_inference() -> None:
    _assert(infer_sequence_type("ACGTNN") == "DNA", "DNA inference failed")
    _assert(infer_sequence_type("MKTLLV") == "PROTEIN", "Protein inference failed")


def _test_program_database_compatibility() -> None:
    BlastParameters(program="blastn", database="nt").normalized()
    BlastParameters(program="blastp", database="nr").normalized()

    try:
        BlastParameters(program="blastn", database="nr").normalized()
    except BlastAnalysisError:
        pass
    else:
        raise AssertionError("Incompatible program/database accepted")


def _test_query_builder() -> None:
    query = build_query(
        "ACGT",
        name="Test",
        sequence_type="DNA",
        source="Unit Test",
        region_start=5,
    )
    _assert(query.length == 4, "Query length incorrect")
    _assert(query.region_start == 5, "Query start incorrect")
    _assert(query.region_end == 8, "Query end incorrect")


def _test_fasta_query_loader() -> None:
    with tempfile.TemporaryDirectory(prefix="synapse_fasta_query_test_") as temp_dir:
        fasta = Path(temp_dir) / "query.fasta"
        fasta.write_text(
            ">Query_A\r\nACGTACGT\r\n>Query_B\r\nGGCCGGCC\r\n",
            encoding="utf-8",
        )

        queries = load_fasta_queries(fasta)
        _assert(len(queries) == 2, "FASTA loader did not recover both records")
        _assert(queries[0].name == "Query_A", "First FASTA header was parsed incorrectly")
        _assert(queries[0].sequence == "ACGTACGT", "First FASTA sequence was parsed incorrectly")
        _assert(queries[1].name == "Query_B", "Second FASTA header was parsed incorrectly")
        _assert(queries[1].source == "FASTA File", "FASTA query source was not recorded")


def _test_xml_parser() -> None:
    query = build_query(
        "ATGCGTACGTACGTACGTAC",
        name="Test_Query",
        sequence_type="DNA",
    )
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="TESTRID",
        parameters=params,
        submitted_at="2026-01-01T00:00:00+00:00",
        completed_at="2026-01-01T00:01:00+00:00",
    )

    _assert(run.hit_count == 1, "XML parser did not recover one hit")
    _assert(run.hits[0].accession == "TEST001", "Accession parsing failed")
    _assert(run.hits[0].best_identity_percent == 95.0, "Identity calculation failed")
    _assert(run.hits[0].best_coverage_percent == 100.0, "Coverage calculation failed")


def _test_hit_ranking() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )

    ranked = rank_hits(run.hits)
    _assert(ranked[0].accession == "TEST001", "Ranking failed")


def _test_evidence_profile() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )

    profile = build_evidence_profile(run)[0]
    _assert(profile.identity_percent == 95.0, "Evidence identity failed")
    _assert(profile.query_coverage_percent == 100.0, "Evidence coverage failed")
    _assert(0 <= profile.evidence_score <= 100, "Evidence score out of range")
    _assert("experimental" not in profile.evidence_band.lower(), "Invalid evidence label")


def _test_alignment() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )
    text = alignment_text(run.hits[0])
    _assert("Query" in text and "Subject" in text, "Alignment rendering failed")


def _test_serialization() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )
    payload = json.loads(run_to_json(run))
    _assert(payload["metadata"]["rid"] == "RID", "JSON metadata failed")
    _assert(payload["hits"][0]["accession"] == "TEST001", "JSON hit failed")


def _test_comparison() -> None:
    query = build_query("ACGT", name="A", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run_a = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="A",
        parameters=params,
        submitted_at="now",
    )

    query_b = build_query("ACGT", name="B", sequence_type="DNA")
    run_b = parse_blast_xml(
        SAMPLE_XML,
        query=query_b,
        rid="B",
        parameters=params,
        submitted_at="now",
    )

    comparison = compare_blast_runs(run_a, run_b)
    _assert(comparison.common_accessions == ["TEST001"], "Comparison failed")


def _test_orf_adapter() -> None:
    orf = {
        "id": "ORF-001",
        "sequence": "ATGCGT",
        "start": 101,
    }
    query = build_orf_query(orf)
    _assert(query.name == "ORF-001", "ORF name extraction failed")
    _assert(query.region_start == 101, "ORF coordinate extraction failed")
    _assert(query.region_end == 106, "ORF end coordinate failed")
    _assert(query.source == "Selected ORF", "ORF source metadata failed")


def _test_export_schema() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )
    payload = run_to_dict(run)
    _assert("evidence_profile" in payload, "Evidence profile missing")
    _assert("database_statistics" in payload, "Database statistics missing")
    _assert("interpretation_notes" in payload, "Interpretation notes missing")



def _test_ncbi_client_configuration() -> None:
    client = NCBIBlastClient(tool="Synapse Test Tool", email="test@example.org")
    _assert(client.tool == "Synapse_Test_Tool", "NCBI tool identifier was not normalized")
    _assert(client.email == "test@example.org", "NCBI email configuration was lost")
    _assert(
        client.min_request_interval >= NCBI_MIN_REQUEST_INTERVAL_SECONDS,
        "NCBI request interval is below the documented minimum",
    )
    _assert(
        client.min_poll_interval >= NCBI_MIN_POLL_INTERVAL_SECONDS,
        "NCBI polling interval is below the documented minimum",
    )


def _test_ncbi_request_parameters() -> None:
    client = NCBIBlastClient(tool="Synapse", email="dev@example.org")
    query = build_query("ACGT", name="NCBI_Query", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    normalized = params.normalized()
    payload = {
        "CMD": "Put",
        "PROGRAM": normalized.program,
        "DATABASE": normalized.database,
        "QUERY": f">{query.name}\r\n{query.sequence}",
        "EXPECT": normalized.evalue,
        "HITLIST_SIZE": normalized.max_target_seqs,
        "TOOL": client.tool,
        "EMAIL": client.email,
    }
    _assert(payload["CMD"] == "Put", "NCBI command parameter missing")
    _assert(payload["TOOL"] == "Synapse", "NCBI tool parameter missing")
    _assert(payload["EMAIL"] == "dev@example.org", "NCBI email parameter missing")
    _assert(payload["DATABASE"] == "nt", "NCBI database parameter incorrect")


def _test_local_database_command_reference() -> None:
    """Ensure local BLAST uses a basename/cwd pair for robust Windows paths."""
    with tempfile.TemporaryDirectory(prefix="Synapse Project With Spaces ") as temp_dir:
        database = Path(temp_dir) / "nested folder" / "synapse_test"
        database.parent.mkdir(parents=True, exist_ok=True)
        reference, working_directory = _build_local_database_command_reference(database)
        _assert(reference == "synapse_test", "Database reference should be the prefix name")
        _assert(
            working_directory == database.parent.resolve(),
            "Database working directory was resolved incorrectly",
        )


def _test_local_blast_environment() -> None:
    """Ensure BLASTDB includes the database directory without losing existing paths."""
    with tempfile.TemporaryDirectory(prefix="synapse_env_test_") as temp_dir:
        database = Path(temp_dir) / "db" / "testdb"
        env = _build_local_blast_environment(database)
        _assert(
            env["BLASTDB"] == str(database.parent.resolve()),
            "BLASTDB was not initialized with the database directory",
        )


def _test_local_blast_environment_preserves_existing() -> None:
    """Ensure an existing BLASTDB search path is preserved."""
    with tempfile.TemporaryDirectory(prefix="synapse_env_test_") as temp_dir:
        database = Path(temp_dir) / "db" / "testdb"
        old_value = os.environ.get("BLASTDB")
        try:
            os.environ["BLASTDB"] = "existing_database_root"
            env = _build_local_blast_environment(database)
        finally:
            if old_value is None:
                os.environ.pop("BLASTDB", None)
            else:
                os.environ["BLASTDB"] = old_value

        parts = env["BLASTDB"].split(os.pathsep)
        _assert(str(database.parent.resolve()) in parts, "Database root was not prepended")
        _assert("existing_database_root" in parts, "Existing BLASTDB path was lost")


def _test_default_database_root_resolution() -> None:
    """Ensure the default root is a project-level data/blast_db path."""
    expected_name = Path("data") / "blast_db"
    resolved = DEFAULT_LOCAL_DB_ROOT.resolve()
    _assert(
        resolved.name == expected_name.name,
        "Default BLAST root must end in data/blast_db",
    )
    _assert(
        resolved.parent.name == expected_name.parent.name,
        "Default BLAST root parent must be data",
    )


def _test_local_database_structure() -> None:
    with tempfile.TemporaryDirectory(prefix="synapse_db_test_") as temp_dir:
        root = Path(temp_dir)
        (root / "testdb.nhr").write_bytes(b"nhr")
        (root / "testdb.nin").write_bytes(b"nin")
        (root / "testdb.nsq").write_bytes(b"nsq")

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("testdb", "nucleotide")

        _assert(database.installed, "Local database was not detected")
        _assert(database.structurally_valid, "Local database structure was not valid")
        _assert(database.file_count == 3, "Unexpected local database file count")


def _test_nested_runtime_alias() -> None:
    """Verify runtime aliases point to sibling nested volume prefixes."""
    with tempfile.TemporaryDirectory(prefix="synapse_runtime_alias_test_") as temp_dir:
        root = Path(temp_dir)
        container = root / "refseq_rna"
        first_volume = container / "refseq_rna.00"
        second_volume = container / "refseq_rna.01"
        first_volume.mkdir(parents=True)
        second_volume.mkdir(parents=True)

        for volume in (first_volume, second_volume):
            prefix = volume / volume.name
            for extension in (".nhr", ".nin", ".nsq"):
                (Path(str(prefix) + extension)).write_bytes(b"test")

        alias = first_volume / "refseq_rna.nal"
        alias.write_text(
            "DBLIST refseq_rna.00 refseq_rna.01\r\n"
            "TITLE Synapse nested runtime alias test\r\n",
            encoding="utf-8",
        )

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("refseq_rna", "nucleotide")
        _assert(database.structurally_valid, "Nested runtime-alias test database was not valid")
        _assert(_needs_runtime_blast_alias(database), "Nested alias did not require runtime alias")

        with tempfile.TemporaryDirectory(prefix="synapse_runtime_alias_out_", dir=str(container)) as runtime_dir:
            reference, working_directory = _create_runtime_blast_alias(
                database, runtime_dir
            )
            runtime_alias = Path(working_directory) / (reference + ".nal")
            payload = runtime_alias.read_text(encoding="utf-8")
            _assert(runtime_alias.suffix == ".nal", "Runtime alias extension was incorrect")
            _assert(reference == "refseq_rna", "Runtime BLAST reference must use the canonical database root")
            _assert("DBLIST refseq_rna.00 refseq_rna.01" in payload, "Runtime alias did not point to volume 00")
            _assert("DBLIST refseq_rna.00 refseq_rna.01" in payload, "Runtime alias did not point to volume 01")
            _assert("DBLIST refseq_rna.00 refseq_rna.01" in payload, "Runtime DBLIST was malformed")
            _assert((Path(working_directory) / "refseq_rna.00.nsq").exists(), "Runtime flat view is missing volume 00")
            _assert((Path(working_directory) / "refseq_rna.01.nsq").exists(), "Runtime flat view is missing volume 01")
            _assert((Path(working_directory) / "refseq_rna.00.nsq").stat().st_ino == (first_volume / "refseq_rna.00.nsq").stat().st_ino, "Runtime volume 00 was copied instead of hard-linked")


def _test_blastdbcmd_info_parser() -> None:
    sample = """Database: Synapse Local Test Database
        3 sequences; 153 total bases

BLASTDB Version: 5

Volumes:
        C:\\Synapse\\synapse_test
"""
    info = _parse_blastdbcmd_info(sample)
    _assert(info["title"] == "Synapse Local Test Database", "Database title parsing failed")
    _assert(info["sequence_count"] == 3, "Sequence count parsing failed")
    _assert(info["total_length"] == 153, "Total length parsing failed")
    _assert(info["blastdb_version"] == "5", "BLASTDB version parsing failed")
    _assert(len(info["volumes"]) == 1, "Volume parsing failed")


def _test_nested_local_database_discovery() -> None:
    """Ensure Synapse detects BLAST DB v5-style files in subdirectories."""
    with tempfile.TemporaryDirectory(prefix="synapse_nested_db_test_") as temp_dir:
        root = Path(temp_dir)
        nested = root / "local_test"
        nested.mkdir(parents=True)

        # BLAST DB v5 commonly exposes a .ndb core file. Additional files make
        # the fixture resemble a real nucleotide database layout.
        (nested / "synapse_test.ndb").write_bytes(b"ndb")
        (nested / "synapse_test.nhr").write_bytes(b"nhr")
        (nested / "synapse_test.nin").write_bytes(b"nin")
        (nested / "synapse_test.nsq").write_bytes(b"nsq")

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("synapse_test")

        _assert(database.installed, "Nested BLAST database was not detected")
        _assert(database.structurally_valid, "BLAST DB v5 structure was not recognized")
        _assert(database.name == "local_test/synapse_test", "Nested database name was not resolved")
        _assert(
            Path(database.path) == (nested / "synapse_test").resolve(),
            "Nested database prefix was resolved incorrectly",
        )

        listed = manager.list_databases()
        _assert(
            any(item.name == "local_test/synapse_test" for item in listed),
            "Nested database was not returned by list_databases",
        )


def _test_nested_multivolume_database_discovery() -> None:
    """Verify the extracted NCBI-style folder-per-volume layout."""
    with tempfile.TemporaryDirectory(prefix="synapse_multivolume_db_test_") as temp_dir:
        root = Path(temp_dir)
        container = root / "refseq_rna"
        container.mkdir(parents=True)
        alias = container / "refseq_rna.nal"
        alias.write_text(
            "DBLIST refseq_rna.00 refseq_rna.01\r\n"
            "TITLE Synapse RefSeq RNA test\r\n",
            encoding="utf-8",
        )

        for volume in ("00", "01"):
            folder = container / f"refseq_rna.{volume}"
            folder.mkdir()
            for extension in (".nhr", ".nin", ".nsq"):
                (folder / f"refseq_rna.{volume}{extension}").write_bytes(b"test")

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("refseq_rna", "nucleotide")

        _assert(database.installed, "Nested multi-volume database was not detected")
        _assert(database.structurally_valid, "Complete multi-volume database was not recognized")
        _assert(database.name == "refseq_rna", "Logical database name changed unexpectedly")
        _assert(database.volume_count == 2, "Expected two BLAST volumes")
        _assert(database.expected_volume_count == 2, "Alias DBLIST volume count was not parsed")
        _assert(not database.missing_volumes, "Complete database reported missing volumes")
        _assert(not database.incomplete_volumes, "Complete volumes were marked incomplete")
        _assert(
            all(Path(path).is_dir() for path in database.volume_paths),
            "Volume directory paths were not discovered correctly",
        )

        # Alias numbering can be unpadded while the extracted filesystem
        # uses the canonical two-digit NCBI volume directories. The resolver
        # must match .0 to .00 (and .1 to .01) by numeric index.
        alias.write_text(
            "DBLIST refseq_rna.0 refseq_rna.1\r\n"
            "TITLE Synapse RefSeq RNA test\r\n",
            encoding="utf-8",
        )
        unpadded = manager.scan_database("refseq_rna", "nucleotide")
        _assert(unpadded.structurally_valid, "Unpadded .0/.1 alias did not match .00/.01 volumes")
        _assert(unpadded.volume_count == 2, "Unpadded alias volume count was incorrect")
        _assert(not unpadded.missing_volumes, "Unpadded alias incorrectly reported missing volumes")
        _assert(not unpadded.incomplete_volumes, "Unpadded alias marked complete volumes incomplete")

        # Restore the canonical padded alias for the remaining damage/missing
        # volume assertions.
        alias.write_text(
            "DBLIST refseq_rna.00 refseq_rna.01\r\n"
            "TITLE Synapse RefSeq RNA test\r\n",
            encoding="utf-8",
        )

        # A misleading prefix inside a volume directory must never replace
        # the parent database alias.
        decoy = container / "refseq_rna.00" / "refseq_rna.nhr"
        decoy.write_bytes(b"decoy")
        resolved = manager.scan_database("refseq_rna", "nucleotide")
        _assert(
            Path(resolved.path) == (container / "refseq_rna").resolve(),
            "A volume-internal prefix incorrectly became the database prefix",
        )

        # Remove one required file to ensure the manager refuses to call the
        # database structurally complete when an extracted volume is damaged.
        (container / "refseq_rna.00" / "refseq_rna.00.nsq").unlink()
        damaged = manager.scan_database("refseq_rna", "nucleotide")
        _assert(not damaged.structurally_valid, "Damaged volume was accepted as valid")
        _assert(
            "refseq_rna.00" in damaged.incomplete_volumes,
            "Damaged volume was not reported",
        )
        _assert(
            damaged.missing_volume_files.get("refseq_rna.00") == [".nsq"],
            "Missing volume file was not reported precisely",
        )

        # Remove the second volume completely to verify missing-volume reporting.
        for path in (container / "refseq_rna.01").iterdir():
            path.unlink()
        (container / "refseq_rna.01").rmdir()
        missing = manager.scan_database("refseq_rna", "nucleotide")
        _assert(
            "refseq_rna.01" in missing.missing_volumes,
            "Missing BLAST volume was not reported",
        )


def _test_alias_inside_first_volume_directory() -> None:
    """Verify an alias stored inside .00 can still expose all sibling volumes."""
    with tempfile.TemporaryDirectory(prefix="synapse_alias_in_first_volume_") as temp_dir:
        root = Path(temp_dir)
        container = root / "refseq_rna"
        container.mkdir(parents=True)
        first = container / "refseq_rna.00"
        first.mkdir()
        (first / "refseq_rna.nal").write_text(
            "DBLIST refseq_rna.00 refseq_rna.01\r\nTITLE Test\r\n",
            encoding="utf-8",
        )
        # Shared alias metadata may also be stored beside the alias.
        (first / "refseq_rna.ndb").write_bytes(b"ndb")
        for volume in ("00", "01"):
            folder = container / f"refseq_rna.{volume}"
            if not folder.exists():
                folder.mkdir()
            for extension in (".nhr", ".nin", ".nsq"):
                (folder / f"refseq_rna.{volume}{extension}").write_bytes(b"test")

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("refseq_rna", "nucleotide")
        _assert(database.structurally_valid, "Alias inside .00 was not recognized as complete")
        _assert(database.volume_count == 2, "Sibling volumes were not discovered")
        _assert(database.expected_volume_count == 2, "Alias DBLIST count was not preserved")
        _assert(not database.missing_volumes, "Sibling volume was incorrectly reported missing")
        _assert(not database.incomplete_volumes, "Complete sibling volumes were marked incomplete")


def _test_local_database_manifest() -> None:
    with tempfile.TemporaryDirectory(prefix="synapse_manifest_test_") as temp_dir:
        root = Path(temp_dir)
        manager = LocalBlastDatabaseManager(root=root)
        manager.ensure_root()

        manifest = {
            "custom_db": {
                "database_type": "protein",
                "source": "Custom FASTA",
            }
        }
        manager._save_manifest(manifest)
        loaded = manager._load_manifest()

        _assert(loaded["custom_db"]["database_type"] == "protein", "Manifest roundtrip failed")


def _test_local_parameter_compatibility() -> None:
    manager = LocalBlastDatabaseManager()
    with tempfile.TemporaryDirectory(prefix="synapse_db_test_") as temp_dir:
        root = Path(temp_dir)
        (root / "protein_db.phr").write_bytes(b"phr")
        (root / "protein_db.pin").write_bytes(b"pin")
        (root / "protein_db.psq").write_bytes(b"psq")
        database = manager.__class__(root=root).scan_database("protein_db", "protein")

        query = build_query("MKTLLV", name="Protein_Query", sequence_type="PROTEIN")
        params = BlastParameters(program="blastp", database="protein_db")

        _assert(
            PROGRAMS[params.program]["database_type"] == database.database_type,
            "Program/database type compatibility failed",
        )
        _assert(query.sequence_type == "PROTEIN", "Protein query type failed")


def _test_synapse_user_path_resolution() -> None:
    application_root = Path(__file__).resolve().parent.parent
    relative_input = Path("relative_test") / "result.json"
    relative = _resolve_synapse_user_path(relative_input)
    system_root = Path(os.environ.get("SystemRoot", r"C:\Windows")).resolve()
    expected_base = application_root if (
        os.name == "nt" and Path.cwd().resolve() == system_root / "System32"
    ) else Path.cwd().resolve()
    _assert(
        relative == (expected_base / relative_input).resolve(),
        "Relative user path was resolved incorrectly",
    )

    absolute = _resolve_synapse_user_path(application_root / "data")
    _assert(absolute == (application_root / "data").resolve(),
            "Absolute path was modified incorrectly")


def _test_executable_manager() -> None:
    manager = BlastExecutableManager()
    available = manager.available()

    _assert(set(BLAST_EXECUTABLES) == set(available), "Executable catalog mismatch")
    _assert(manager.resolve_perl() is None or isinstance(manager.resolve_perl(), str),
            "Perl discovery returned an invalid value")


def _test_short_query_rescue_policy() -> None:
    query = build_query(
        "ATGAAACGTGCTGACTTTGGTGACGCTTACCGTAAAGGTC" * 2,
        name="Short_DNA_Query",
        sequence_type="DNA",
    )
    default_parameters = BlastParameters(
        program="blastn",
        database="nt",
        evalue=10.0,
    )
    _assert(
        _should_run_short_query_rescue(query, default_parameters),
        "Default short-query rescue policy should be eligible",
    )

    explicit_word_size = BlastParameters(
        program="blastn",
        database="nt",
        word_size=7,
    )
    _assert(
        not _should_run_short_query_rescue(query, explicit_word_size),
        "Explicit word size must disable automatic rescue",
    )


def _test_short_query_task_configuration() -> None:
    parameters = BlastParameters(
        program="blastn",
        database="nt",
        task="blastn-short",
        filter_query="no",
    ).normalized()
    _assert(parameters.task == "blastn-short", "Short-query task was not preserved")
    _assert(parameters.filter_query == "no", "Short-query filtering mode was not preserved")


def _test_taxonomy_definition_fallback() -> None:
    name = _infer_scientific_name_from_definition(
        "PREDICTED: Nicrophorus vespilloides protein UBASH3A homolog"
    )
    _assert(name == "Nicrophorus vespilloides", "Failed to infer binomial organism name")

    bracketed = _infer_scientific_name_from_definition(
        "Reference transcript [Homo sapiens]"
    )
    _assert(bracketed == "Homo sapiens", "Failed to parse bracketed organism name")

    _assert(
        _infer_scientific_name_from_definition("PREDICTED: hypothetical protein") is None,
        "Should not infer a non-taxonomic phrase",
    )


def _test_taxonomy_grouping_with_fallback() -> None:
    hit_a = BlastHit(
        hit_id="A",
        accession="A",
        definition="PREDICTED: Nicrophorus vespilloides protein UBASH3A homolog",
        length=100,
    )
    hit_b = BlastHit(
        hit_id="B",
        accession="B",
        definition="PREDICTED: Nicrophorus vespilloides transcript variant",
        length=100,
    )
    hit_c = BlastHit(
        hit_id="C",
        accession="C",
        definition="Example reference sequence",
        length=100,
    )
    for hit in (hit_a, hit_b, hit_c):
        inferred = _infer_scientific_name_from_definition(hit.definition)
        if inferred:
            hit.scientific_names.append(inferred)

    groups = group_hits_by_taxonomy((hit_a, hit_b, hit_c))
    _assert(groups.get("Nicrophorus vespilloides") == 2, "Taxonomy grouping failed")
    _assert(groups.get("Unknown / not provided") == 1, "Unknown taxonomy count failed")


def run_self_test(verbose: bool = True) -> bool:
    tests = [
        ("sequence_validation", _test_sequence_validation),
        ("type_inference", _test_type_inference),
        ("program_database_compatibility", _test_program_database_compatibility),
        ("query_builder", _test_query_builder),
        ("fasta_query_loader", _test_fasta_query_loader),
        ("xml_parser", _test_xml_parser),
        ("hit_ranking", _test_hit_ranking),
        ("evidence_profile", _test_evidence_profile),
        ("alignment_renderer", _test_alignment),
        ("json_serialization", _test_serialization),
        ("run_comparison", _test_comparison),
        ("orf_adapter", _test_orf_adapter),
        ("export_schema", _test_export_schema),
        ("ncbi_client_configuration", _test_ncbi_client_configuration),
        ("ncbi_request_parameters", _test_ncbi_request_parameters),
        ("local_database_command_reference", _test_local_database_command_reference),
        ("local_blast_environment", _test_local_blast_environment),
        ("local_blast_environment_preserves_existing", _test_local_blast_environment_preserves_existing),
        ("default_database_root_resolution", _test_default_database_root_resolution),
        ("local_database_structure", _test_local_database_structure),
        ("nested_local_database_discovery", _test_nested_local_database_discovery),
        ("nested_multivolume_database_discovery", _test_nested_multivolume_database_discovery),
        ("nested_runtime_alias", _test_nested_runtime_alias),
        ("alias_inside_first_volume_directory", _test_alias_inside_first_volume_directory),
        ("blastdbcmd_info_parser", _test_blastdbcmd_info_parser),
        ("local_database_manifest", _test_local_database_manifest),
        ("local_parameter_compatibility", _test_local_parameter_compatibility),
        ("synapse_user_path_resolution", _test_synapse_user_path_resolution),
        ("executable_manager", _test_executable_manager),
        ("short_query_rescue_policy", _test_short_query_rescue_policy),
        ("short_query_task_configuration", _test_short_query_task_configuration),
        ("taxonomy_definition_fallback", _test_taxonomy_definition_fallback),
        ("taxonomy_grouping_with_fallback", _test_taxonomy_grouping_with_fallback),
    ]

    results: List[SelfTestResult] = []

    for name, test in tests:
        try:
            test()
            results.append(SelfTestResult(name, True))
        except Exception as exc:
            results.append(SelfTestResult(name, False, str(exc)))

    if verbose:
        print("\r\n" + _line("="))
        print("SYNAPSE BLAST SELF-TEST")
        print(_line("="))
        for result in results:
            marker = "PASS" if result.passed else "FAIL"
            suffix = f" - {result.detail}" if result.detail else ""
            print(f"[{marker}] {result.name}{suffix}")
        passed = sum(result.passed for result in results)
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
    print("Synapse BLAST module is ready.")
    print("No network search was executed during self-test.")
    print("Use blast_menu(...) from main.py for interactive BLAST Search.")
