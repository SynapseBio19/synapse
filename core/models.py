"""Shared data model: sequences, intervals, features and analysis results.

Design rules (audit Phase III/V):

* One internal coordinate convention: 0-based, half-open ``[start, end)`` on
  the FORWARD strand of the parent sequence, plus an explicit strand.
  Human-facing output uses ``to_1based()`` (1-based, inclusive).
* Every stored result names the module, its parameters and the input it was
  computed from (sequence id + SHA-256), so a report can be reproduced.
* The model is additive: existing modules keep returning their dictionaries;
  ``AnalysisResult.data`` holds them unchanged.
"""

from __future__ import annotations

import hashlib
import platform
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from core.alphabets import classify_sequence

SYNAPSE_VERSION = "0.2.0.dev0"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Interval:
    """0-based, half-open interval on the forward strand, with a strand."""

    start: int
    end: int
    strand: str = "+"

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < self.start:
            raise ValueError(f"Invalid interval [{self.start}, {self.end}).")
        if self.strand not in {"+", "-", "."}:
            raise ValueError("strand must be '+', '-' or '.'.")

    @classmethod
    def from_1based(cls, start: int, end: int, strand: str = "+") -> "Interval":
        """From 1-based inclusive coordinates (as shown to users / in VCF)."""
        return cls(start - 1, end, strand)

    @classmethod
    def from_0based_inclusive(cls, start: int, stop: int, strand: str = "+") -> "Interval":
        """From 0-based inclusive coordinates (Synapse ORF dict convention)."""
        return cls(start, stop + 1, strand)

    def to_1based(self) -> Tuple[int, int]:
        return self.start + 1, self.end

    @property
    def length(self) -> int:
        return self.end - self.start

    def extract(self, parent: str) -> str:
        """Sequence of the interval in its own 5'->3' orientation."""
        from core.seqtools import reverse_complement

        if self.end > len(parent):
            raise ValueError("Interval extends beyond the parent sequence.")
        segment = parent[self.start:self.end]
        return reverse_complement(segment) if self.strand == "-" else segment

    def overlaps(self, other: "Interval") -> bool:
        return self.start < other.end and other.start < self.end

    def describe(self) -> str:
        start, end = self.to_1based()
        return f"{start}-{end} ({self.strand})"


# ---------------------------------------------------------------------------
# Sequence record and features
# ---------------------------------------------------------------------------

@dataclass
class Feature:
    kind: str
    interval: Interval
    name: str = ""
    qualifiers: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SequenceRecord:
    """A loaded sequence with its provenance and sequence type.

    ``alphabet`` (DNA, RNA, PROTEIN, or None when undetermined) and
    ``type_evidence`` come from ``core.alphabets.classify_sequence``.
    ``declared_type`` is a type stated by the user or by file metadata;
    an ``alphabet`` passed to the constructor counts as a declaration.
    ``molecule_type`` is the file's own label (e.g. GenBank "mRNA").
    """

    identifier: str
    sequence: str
    description: str = ""
    topology: str = "linear"
    alphabet: Optional[str] = None
    source_path: Optional[str] = None
    source_format: Optional[str] = None
    source_sha256: Optional[str] = None
    record_index: Optional[int] = None
    features: List[Feature] = field(default_factory=list)
    declared_type: Optional[str] = None
    molecule_type: Optional[str] = None
    type_evidence: str = field(default="", init=False)
    type_reasons: Tuple[str, ...] = field(default=(), init=False)

    def __post_init__(self) -> None:
        self.sequence = "".join(str(self.sequence).split()).upper()
        if self.topology not in {"linear", "circular"}:
            raise ValueError("topology must be 'linear' or 'circular'.")
        if self.declared_type is None and self.alphabet is not None:
            self.declared_type = self.alphabet
        result = classify_sequence(self.sequence, self.declared_type)
        self.declared_type = result.declared
        self.alphabet = result.sequence_type
        self.type_evidence = result.evidence
        self.type_reasons = result.reasons

    def with_declared_type(self, declared_type: str) -> "SequenceRecord":
        """A copy with *declared_type* (DNA, RNA, PROTEIN) and its classification."""
        return replace(self, alphabet=None, declared_type=declared_type)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sequence_sha256(self) -> str:
        return sha256_text(self.sequence)

    def provenance(self) -> Dict[str, Any]:
        return {
            "identifier": self.identifier,
            "description": self.description,
            "length": self.length,
            "topology": self.topology,
            "sequence_type": self.alphabet,
            "sequence_type_evidence": self.type_evidence,
            "declared_type": self.declared_type,
            "molecule_type": self.molecule_type,
            "sequence_sha256": self.sequence_sha256,
            "source_path": self.source_path,
            "source_format": self.source_format,
            "source_sha256": self.source_sha256,
            "record_index": self.record_index,
        }

    @classmethod
    def from_legacy_dict(cls, record: Dict[str, Any], index: Optional[int] = None) -> "SequenceRecord":
        """Build from the dicts returned by ``fileio`` loaders."""
        header = str(record.get("header", "")).strip()
        identifier = header.split()[0] if header else f"record_{index or 1}"
        return cls(
            identifier=identifier,
            sequence=record.get("sequence", ""),
            description=header,
            topology=record.get("topology", "linear"),
            declared_type=record.get("declared_type"),
            molecule_type=record.get("molecule_type"),
            source_path=record.get("source_path"),
            source_format=record.get("source_format"),
            source_sha256=record.get("source_sha256"),
            record_index=index,
        )


# ---------------------------------------------------------------------------
# Analysis results
# ---------------------------------------------------------------------------

@dataclass
class AnalysisResult:
    """One stored analysis with everything needed to reproduce it."""

    module: str
    data: Any
    parameters: Dict[str, Any] = field(default_factory=dict)
    input_sequence_sha256: Optional[str] = None
    input_description: str = ""
    warnings: List[str] = field(default_factory=list)
    limitations: List[str] = field(default_factory=list)
    created_at: str = field(default_factory=utc_now)
    synapse_version: str = SYNAPSE_VERSION
    python_version: str = field(default_factory=platform.python_version)

    def summary(self) -> Dict[str, Any]:
        """Metadata without the (possibly large) data payload."""
        payload = asdict(self)
        payload.pop("data", None)
        return payload
