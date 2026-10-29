"""
Synapse - Variant Analysis
==========================

Professional variant-level analysis layer for Project Synapse.

Scope
-----
This module is deliberately separated from mutation.py.

mutation.py:
    sequence-level mutation detection and molecular consequence analysis.

variant_analysis.py:
    variant records, VCF-like input, filtering, normalization, annotation
    containers, cohort analysis, variant comparison, frequencies, quality
    metrics, summaries, reporting, JSON/TSV export, and CLI.

The module is self-contained and does not require third-party packages.
When mutation.py is available, sequence-derived consequence analysis can be
connected through the optional adapter functions near the end of the file.

Important scientific boundary
-----------------------------
This module does not pretend to provide clinical interpretation or database
annotation without evidence. Fields such as clinical significance,
population frequency, gene, transcript, and pathogenicity are represented
explicitly and are only populated when supplied by the user or an external
annotation source.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
import csv
import io
import json
import math
import re
import statistics
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


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


# VCF 4.3 reserved field cardinalities (spec sections 1.6.1 and 1.6.2). Used
# when the file header does not declare a field.
RESERVED_INFO_NUMBER = {
    "AA": "1", "AC": "A", "AD": "R", "ADF": "R", "ADR": "R", "AF": "A",
    "AN": "1", "BQ": "1", "CIGAR": "A", "DB": "0", "DP": "1", "END": "1",
    "H2": "0", "H3": "0", "MQ": "1", "MQ0": "1", "NS": "1", "SB": "4",
    "SOMATIC": "0", "VALIDATED": "0", "1000G": "0",
}
RESERVED_FORMAT_NUMBER = {
    "AD": "R", "ADF": "R", "ADR": "R", "DP": "1", "EC": "A", "FT": "1",
    "GL": "G", "GP": "G", "GQ": "1", "GT": "1", "HQ": "2", "MQ": "1",
    "PL": "G", "PQ": "1", "PS": "1",
}

_HEADER_FIELD_RE = re.compile(r"^##(INFO|FORMAT)=<(.*)>$")


def parse_vcf_header_numbers(metadata: Sequence[str]) -> Dict[str, Dict[str, str]]:
    """Return {"INFO": {ID: Number}, "FORMAT": {ID: Number}} from ## header lines."""
    numbers: Dict[str, Dict[str, str]] = {"INFO": {}, "FORMAT": {}}
    for line in metadata:
        match = _HEADER_FIELD_RE.match(line.strip())
        if not match:
            continue
        section, body = match.groups()
        fields = dict(
            (key.strip(), value.strip())
            for key, _, value in (
                part.partition("=") for part in re.split(r",(?=[A-Za-z]+=)", body)
            )
        )
        if "ID" in fields and "Number" in fields:
            numbers[section][fields["ID"]] = fields["Number"]
    return numbers


def _field_number(section: str, key: str, header_numbers: Optional[Mapping[str, Mapping[str, str]]]) -> Optional[str]:
    if header_numbers and key in header_numbers.get(section, {}):
        return header_numbers[section][key]
    reserved = RESERVED_INFO_NUMBER if section == "INFO" else RESERVED_FORMAT_NUMBER
    return reserved.get(key)


def _genotype_index(j: int, k: int) -> int:
    """VCF ordering of diploid genotype j/k (j <= k): F(j,k) = k(k+1)/2 + j."""
    return k * (k + 1) // 2 + j


def _select_allele_specific_info(
    info: Mapping[str, Any],
    alt_index: int,
    alt_count: int,
    header_numbers: Optional[Mapping[str, Mapping[str, str]]] = None,
) -> Dict[str, Any]:
    """Select allele-specific INFO values when splitting a multi-allelic record.

    Cardinality comes from the ``##INFO`` header or, for reserved keys, the
    VCF 4.3 specification: Number=A keeps the selected ALT's value, Number=R
    keeps [REF, selected ALT]. Undeclared list fields fall back to the
    previous heuristic (list length == ALT count -> treated as Number=A).
    """
    selected: Dict[str, Any] = {}
    for key, value in info.items():
        if not isinstance(value, (list, tuple)):
            selected[key] = value
            continue
        values = list(value)
        number = _field_number("INFO", key, header_numbers)
        if number == "A" and len(values) == alt_count:
            selected[key] = values[alt_index - 1]
        elif number == "R" and len(values) == alt_count + 1:
            selected[key] = [values[0], values[alt_index]]
        elif number is None and len(values) == alt_count:
            selected[key] = values[alt_index - 1]
        else:
            selected[key] = values
    return selected


def _remap_sample_data_for_alt(
    values: Mapping[str, str],
    alt_index: int,
    alt_count: Optional[int] = None,
    header_numbers: Optional[Mapping[str, Mapping[str, str]]] = None,
) -> Dict[str, str]:
    """Remap GT and allele-specific sample fields for one ALT allele.

    GT: the selected ALT becomes 1; REF stays 0; every other ALT is shown as
    0 (REF-masked, the same convention as ``bcftools norm -m-``).
    Number=R fields (e.g. AD) keep [REF, ALT]; Number=A keep the ALT value;
    Number=G fields (e.g. PL, GL, GP) keep the genotypes 0/0, 0/a, a/a for a
    diploid call (or 0, a for a haploid call), using the VCF genotype order.
    """
    remapped = dict(values)
    ploidy = None
    gt = values.get("GT")
    if gt not in (None, "", "."):
        tokens = GT_ALLELE_SEPARATOR_RE.split(gt)
        ploidy = len(tokens)
        separator = "|" if "|" in gt else "/"
        new_tokens: List[str] = []
        for token in tokens:
            if token in ("", "."):
                new_tokens.append(".")
            else:
                try:
                    allele = int(token)
                except ValueError:
                    new_tokens.append(".")
                    continue
                if allele == 0:
                    new_tokens.append("0")
                elif allele == alt_index:
                    new_tokens.append("1")
                else:
                    new_tokens.append("0")
        remapped["GT"] = separator.join(new_tokens)

    for key, raw in values.items():
        if key == "GT" or not isinstance(raw, str) or raw in ("", "."):
            continue
        number = _field_number("FORMAT", key, header_numbers)
        parts = raw.split(",")
        if number == "R" and len(parts) > alt_index:
            remapped[key] = f"{parts[0]},{parts[alt_index]}"
        elif number == "A" and alt_count and len(parts) == alt_count:
            remapped[key] = parts[alt_index - 1]
        elif number == "G" and alt_count:
            if ploidy == 1 and len(parts) == alt_count + 1:
                remapped[key] = f"{parts[0]},{parts[alt_index]}"
            elif ploidy in (None, 2) and len(parts) == _genotype_index(alt_count, alt_count) + 1:
                keep = (
                    _genotype_index(0, 0),
                    _genotype_index(0, alt_index),
                    _genotype_index(alt_index, alt_index),
                )
                remapped[key] = ",".join(parts[i] for i in keep)

    return remapped


@dataclass
class VariantRecord:
    """A VCF-style multi-allelic record."""

    chromosome: str
    position: int
    reference: str
    alternatives: List[str]
    variant_id: Optional[str] = None
    quality: Optional[float] = None
    filter_status: str = "UNFILTERED"
    info: Dict[str, Any] = field(default_factory=dict)
    format_keys: List[str] = field(default_factory=list)
    sample_data: Dict[str, Dict[str, str]] = field(default_factory=dict)
    header_numbers: Optional[Dict[str, Dict[str, str]]] = None

    def split(self) -> List[Variant]:
        """Split a multi-allelic record into independent biallelic variants.

        Allele-specific INFO values (for example ``AF=A,B``) are reduced to
        the value belonging to the selected ALT allele. Genotypes are remapped
        to the new biallelic allele space so that cohort analysis remains
        correct after splitting.
        """
        variants: List[Variant] = []
        for alt_index, alt in enumerate(self.alternatives, start=1):
            if not alt or alt == ".":
                continue
            info = _select_allele_specific_info(
                self.info, alt_index, len(self.alternatives), self.header_numbers
            )
            sample_data = {
                sample: _remap_sample_data_for_alt(
                    values, alt_index, len(self.alternatives), self.header_numbers
                )
                for sample, values in self.sample_data.items()
            }
            variants.append(
                Variant(
                    chromosome=self.chromosome,
                    position=self.position,
                    reference=self.reference,
                    alternative=alt,
                    variant_id=self.variant_id,
                    quality=self.quality,
                    filter_status=self.filter_status,
                    info=info,
                    format_keys=list(self.format_keys),
                    sample_data=sample_data,
                )
            )
        return variants


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


# ============================================================================
# Parsing / validation
# ============================================================================

def parse_float(value: Any) -> Optional[float]:
    if value in (None, "", ".", "NA", "N/A"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        raise VariantAnalysisError(f"Invalid numeric value: {value!r}")


def parse_int(value: Any) -> Optional[int]:
    if value in (None, "", ".", "NA", "N/A"):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        raise VariantAnalysisError(f"Invalid integer value: {value!r}")


def validate_allele(allele: str, *, allow_symbolic: bool = True) -> bool:
    allele = str(allele).upper()
    if allow_symbolic and allele.startswith("<") and allele.endswith(">"):
        return True
    if allele == "*":
        return True
    if not allele:
        return False
    return all(base in DNA_BASES for base in allele)


def validate_variant(variant: Variant) -> List[str]:
    """Return validation warnings/errors without silently modifying a record."""
    issues = []

    if variant.position < 1:
        issues.append("Position must be >= 1.")
    if not validate_allele(variant.reference):
        issues.append("Invalid reference allele.")
    if not validate_allele(variant.alternative):
        issues.append("Invalid alternative allele.")
    if variant.filter_status not in VALID_FILTER_STATUS:
        issues.append(
            f"Unknown FILTER status: {variant.filter_status}."
        )
    if variant.quality is not None and variant.quality < 0:
        issues.append("Quality score is negative.")
    if variant.allele_frequency is not None and not (
        0.0 <= variant.allele_frequency <= 1.0
    ):
        issues.append("Allele frequency must be between 0 and 1.")
    if variant.allele_count is not None and variant.allele_count < 0:
        issues.append("Allele count cannot be negative.")
    if variant.allele_number is not None and variant.allele_number < 0:
        issues.append("Allele number cannot be negative.")
    return issues


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


def parse_info_field(info_text: str) -> Dict[str, Any]:
    """Parse a VCF INFO field into a dictionary."""
    if not info_text or info_text == ".":
        return {}

    info: Dict[str, Any] = {}
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" not in item:
            info[item] = True
            continue
        key, value = item.split("=", 1)
        if "," in value:
            values = value.split(",")
            info[key] = [
                _coerce_info_value(v) for v in values
            ]
        else:
            info[key] = _coerce_info_value(value)
    return info


def _coerce_info_value(value: str) -> Any:
    if value in ("", "."):
        return None
    try:
        if re.fullmatch(r"-?\d+", value):
            return int(value)
        if re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?", value):
            return float(value)
    except ValueError:
        pass
    return value


def format_info_field(info: Mapping[str, Any]) -> str:
    """Serialize an INFO dictionary to VCF-compatible text."""
    if not info:
        return "."

    fields = []
    for key in sorted(info):
        value = info[key]
        if value is True:
            fields.append(str(key))
        elif isinstance(value, (list, tuple)):
            fields.append(f"{key}=" + ",".join(
                "." if item is None else str(item)
                for item in value
            ))
        elif value is None:
            fields.append(f"{key}=.")
        else:
            fields.append(f"{key}={value}")
    return ";".join(fields)


def parse_genotype(
    genotype: str,
) -> Tuple[Optional[List[Optional[int]]], bool]:
    """Parse a GT value and return allele indexes plus phased state."""
    if genotype in (None, "", "."):
        return None, False

    phased = "|" in genotype
    tokens = GT_ALLELE_SEPARATOR_RE.split(genotype)

    alleles: List[Optional[int]] = []
    for token in tokens:
        if token in ("", "."):
            alleles.append(None)
        else:
            try:
                alleles.append(int(token))
            except ValueError:
                return None, phased
    return alleles, phased


def genotype_category(
    genotype: Optional[str],
) -> str:
    """Classify a diploid genotype call."""
    alleles, _ = parse_genotype(genotype)
    if not alleles:
        return "missing"
    if all(a is None for a in alleles):
        return "missing"
    if all(a == 0 for a in alleles if a is not None) and all(
        a is not None for a in alleles
    ):
        return "homozygous_ref"
    if any(a is None for a in alleles):
        return "partial"
    if all(a == alleles[0] for a in alleles) and alleles[0] not in (0, None):
        return "homozygous_alt"
    if any(a not in (0, None) for a in alleles):
        return "heterozygous"
    return "missing"


def _normalize_filter_status(value: Any) -> str:
    """Normalize VCF FILTER notation to Synapse's internal status."""
    if value in (None, "", "."):
        return "UNFILTERED"
    text = str(value).upper()
    if text == "PASS":
        return "PASS"
    if text in {"FAIL", "UNFILTERED", "UNKNOWN"}:
        return text
    return text


def parse_vcf(
    text: str,
    *,
    split_multiallelic: bool = True,
) -> Tuple[List[Variant], List[str], List[str]]:
    """
    Parse VCF text.

    Returns:
        variants, metadata_lines, sample_names
    """
    variants: List[Variant] = []
    metadata: List[str] = []
    sample_names: List[str] = []
    header_numbers: Optional[Dict[str, Dict[str, str]]] = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n\r")
        if not line:
            continue

        if line.startswith("##"):
            metadata.append(line)
            continue

        if line.startswith("#CHROM"):
            columns = line.split("\t")
            if len(columns) > 9:
                sample_names = columns[9:]
            header_numbers = parse_vcf_header_numbers(metadata)
            continue

        if line.startswith("#"):
            continue

        columns = line.split("\t")
        if len(columns) < 8:
            raise VariantAnalysisError(
                "VCF record must contain at least 8 tab-delimited columns."
            )

        chrom, pos, variant_id, ref, alt_text, qual, filt, info_text = columns[:8]
        position = parse_int(pos)
        if position is None:
            raise VariantAnalysisError(f"Invalid VCF position: {pos!r}")

        quality = parse_float(qual)
        alternatives = alt_text.split(",")

        info = parse_info_field(info_text)

        format_keys: List[str] = []
        sample_data: Dict[str, Dict[str, str]] = {}

        if len(columns) > 8:
            format_keys = columns[8].split(":")
            for sample_name, sample_text in zip(sample_names, columns[9:]):
                values = sample_text.split(":")
                sample_data[sample_name] = {
                    key: values[i] if i < len(values) else "."
                    for i, key in enumerate(format_keys)
                }

        record = VariantRecord(
            chromosome=chrom,
            position=position,
            reference=ref,
            alternatives=alternatives,
            variant_id=None if variant_id == "." else variant_id,
            quality=quality,
            filter_status=_normalize_filter_status(filt),
            info=info,
            format_keys=format_keys,
            sample_data=sample_data,
            header_numbers=header_numbers,
        )

        if split_multiallelic:
            variants.extend(record.split())
        else:
            variants.append(
                Variant(
                    chromosome=chrom,
                    position=position,
                    reference=ref,
                    alternative=alt_text,
                    variant_id=record.variant_id,
                    quality=quality,
                    filter_status=_normalize_filter_status(filt),
                    info=info,
                    format_keys=format_keys,
                    sample_data=sample_data,
                )
            )

    return variants, metadata, sample_names


def parse_vcf_file(path: str, *, split_multiallelic: bool = True):
    """Parse a VCF file; ``.gz``/BGZF-compressed files are read transparently."""
    import gzip

    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return parse_vcf(handle.read(), split_multiallelic=split_multiallelic)
    with open(path, "r", encoding="utf-8") as handle:
        return parse_vcf(handle.read(), split_multiallelic=split_multiallelic)


# ============================================================================
# Normalization
# ============================================================================

def trim_common_bases(
    position: int,
    reference: str,
    alternative: str,
) -> Tuple[int, str, str]:
    """
    Trim common suffix/prefix bases while preserving valid VCF-style alleles.

    The function does not left-shift through repeats; use
    left_normalize_variant() for that.
    """
    pos = position
    ref = reference.upper()
    alt = alternative.upper()

    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    return pos, ref, alt


def _reference_allele_matches(allele: str, observed: str) -> bool:
    """REF allele vs. reference bases; N in either matches any base."""
    return all(a == b or a == "N" or b == "N" for a, b in zip(allele, observed))


def left_normalize_variant(
    variant: Variant,
    reference_sequence: Optional[str] = None,
) -> Variant:
    """
    Left-normalize a simple variant.

    A reference sequence is required for repeat-aware left shifting.
    Without it, common-base trimming is still performed.
    """
    new = Variant(
        chromosome=variant.chromosome,
        position=variant.position,
        reference=variant.reference,
        alternative=variant.alternative,
        variant_id=variant.variant_id,
        quality=variant.quality,
        filter_status=variant.filter_status,
        info=dict(variant.info),
        format_keys=list(variant.format_keys),
        sample_data={
            sample: dict(values)
            for sample, values in variant.sample_data.items()
        },
        variant_type=variant.variant_type,
        gene=variant.gene,
        transcript=variant.transcript,
        region=variant.region,
        consequence=variant.consequence,
        protein_change=variant.protein_change,
        coding_change=variant.coding_change,
        allele_frequency=variant.allele_frequency,
        allele_count=variant.allele_count,
        allele_number=variant.allele_number,
        normalized=variant.normalized,
        normalization_notes=list(variant.normalization_notes),
        annotation=dict(variant.annotation),
    )

    if new.variant_type == "SYMBOLIC" or new.alternative == "*" or not validate_allele(
        new.reference, allow_symbolic=False
    ) or not validate_allele(new.alternative, allow_symbolic=False):
        new.normalization_notes.append(
            "Symbolic, spanning-deletion or non-ACGTN allele: not normalized."
        )
        return new

    pos, ref, alt = new.position, new.reference, new.alternative

    reference_checked = False
    if reference_sequence is not None:
        reference_sequence = str(reference_sequence).upper()
        observed = reference_sequence[pos - 1:pos - 1 + len(ref)]
        if len(observed) != len(ref) or not _reference_allele_matches(ref, observed):
            flags = new.annotation.setdefault("validation_flags", [])
            if "REF_MISMATCH" not in flags:
                flags.append("REF_MISMATCH")
            new.normalization_notes.append(
                f"REF allele '{ref}' does not match the reference sequence "
                f"('{observed or 'out of range'}') at position {pos}; "
                "repeat-aware shifting was not applied."
            )
            reference_sequence = None
        else:
            reference_checked = True

    if reference_sequence is not None and len(ref) != len(alt):
        # Tan, Abecasis & Kang (2015), Bioinformatics 31:2202, Algorithm 1:
        # repeat: (1) drop an identical last base from both alleles;
        #         (2) if an allele became empty, prepend the preceding
        #             reference base to both and move one position left;
        # until neither applies. Then trim identical leading bases while
        # both alleles are longer than one base.
        changed = True
        while changed:
            changed = False
            if ref and alt and ref[-1] == alt[-1]:
                ref, alt = ref[:-1], alt[:-1]
                changed = True
            if (not ref or not alt) and pos > 1:
                previous = reference_sequence[pos - 2]
                ref, alt = previous + ref, previous + alt
                pos -= 1
                changed = True
        if (not ref or not alt) and pos == 1:
            # VCF 4.3: an indel at position 1 is anchored on the base AFTER
            # the event.
            following = reference_sequence[len(ref):len(ref) + 1]
            ref, alt = ref + following, alt + following
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref, alt = ref[1:], alt[1:]
            pos += 1
    else:
        pos, ref, alt = trim_common_bases(pos, ref, alt)

    new.position = pos
    new.reference = ref
    new.alternative = alt
    new.variant_type = classify_variant_type(ref, alt)
    new.end_position = pos + max(0, len(ref) - 1)
    new.normalized = True
    new.normalization_notes.append("Common bases trimmed.")
    if reference_checked:
        new.normalization_notes.append(
            "REF allele verified against the reference; left-aligned and "
            "parsimonious (Tan et al. 2015)."
        )
    return new


def normalize_variants(
    variants: Sequence[Variant],
    *,
    reference_sequence: Optional[str] = None,
    contig_references: Optional[Mapping[str, str]] = None,
) -> List[Variant]:
    """Normalize variants.

    ``contig_references`` maps CHROM names to reference sequences; variants on
    an unmapped contig are only trimmed (no repeat-aware shift, no REF check).
    ``reference_sequence`` applies one reference to every variant (legacy use).
    """
    normalized = []
    for v in variants:
        reference = reference_sequence
        if contig_references is not None:
            reference = contig_references.get(v.chromosome)
        normalized.append(left_normalize_variant(v, reference_sequence=reference))
    return normalized


def check_reference_alleles(
    variants: Sequence[Variant],
    contig_references: Mapping[str, str],
) -> Dict[str, Any]:
    """Compare REF alleles with reference sequences (1-based VCF positions).

    Returns counts and the list of mismatching variant keys; mismatches are
    flagged in ``variant.annotation["validation_flags"]``.
    """
    checked = mismatched = unmapped = 0
    mismatches: List[Tuple[str, int, str, str]] = []
    for variant in variants:
        reference = contig_references.get(variant.chromosome)
        if reference is None or not validate_allele(variant.reference, allow_symbolic=False):
            unmapped += 1
            continue
        reference = reference.upper()
        observed = reference[variant.position - 1:variant.position - 1 + len(variant.reference)]
        checked += 1
        if len(observed) != len(variant.reference) or not _reference_allele_matches(variant.reference, observed):
            mismatched += 1
            mismatches.append(variant.key)
            flags = variant.annotation.setdefault("validation_flags", [])
            if "REF_MISMATCH" not in flags:
                flags.append("REF_MISMATCH")
    return {"checked": checked, "mismatched": mismatched, "not_checked": unmapped,
            "mismatches": mismatches}


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
# Filtering
# ============================================================================

def variant_passes_filter(
    variant: Variant,
    criteria: VariantFilter,
) -> Tuple[bool, List[str]]:
    """Evaluate one variant against all configured criteria."""
    reasons: List[str] = []

    if criteria.minimum_quality is not None:
        if variant.quality is None:
            if not criteria.include_unknown_quality:
                reasons.append("Missing quality.")
        elif variant.quality < criteria.minimum_quality:
            reasons.append("Quality below minimum.")

    if criteria.maximum_quality is not None:
        if variant.quality is not None and variant.quality > criteria.maximum_quality:
            reasons.append("Quality above maximum.")

    if criteria.require_pass and variant.filter_status != "PASS":
        reasons.append("FILTER is not PASS.")

    if criteria.allowed_filter_statuses is not None:
        allowed = {s.upper() for s in criteria.allowed_filter_statuses}
        if variant.filter_status.upper() not in allowed:
            reasons.append("FILTER status not allowed.")

    if criteria.variant_types is not None:
        allowed = {s.upper() for s in criteria.variant_types}
        if variant.variant_type.upper() not in allowed:
            reasons.append("Variant type not allowed.")

    if criteria.chromosomes is not None:
        if variant.chromosome not in criteria.chromosomes:
            reasons.append("Chromosome not allowed.")

    if criteria.genes is not None:
        if variant.gene is None:
            reasons.append("Gene annotation missing.")
        elif variant.gene not in criteria.genes:
            reasons.append("Gene not allowed.")

    if criteria.consequences is not None:
        if variant.consequence is None:
            reasons.append("Consequence annotation missing.")
        elif variant.consequence not in criteria.consequences:
            reasons.append("Consequence not allowed.")

    if criteria.minimum_allele_frequency is not None:
        if variant.allele_frequency is None:
            if not criteria.include_unknown_frequency:
                reasons.append("Missing allele frequency.")
        elif variant.allele_frequency < criteria.minimum_allele_frequency:
            reasons.append("Allele frequency below minimum.")

    if criteria.maximum_allele_frequency is not None:
        if variant.allele_frequency is None:
            if not criteria.include_unknown_frequency:
                reasons.append("Missing allele frequency.")
        elif variant.allele_frequency > criteria.maximum_allele_frequency:
            reasons.append("Allele frequency above maximum.")

    if criteria.minimum_allele_count is not None:
        if variant.allele_count is None:
            reasons.append("Missing allele count.")
        elif variant.allele_count < criteria.minimum_allele_count:
            reasons.append("Allele count below minimum.")

    if criteria.maximum_allele_count is not None:
        if variant.allele_count is not None and variant.allele_count > criteria.maximum_allele_count:
            reasons.append("Allele count above maximum.")

    return not reasons, reasons


def filter_variants(
    variants: Sequence[Variant],
    criteria: VariantFilter,
) -> VariantFilterResult:
    variants = synchronize_variant_fields(variants)
    passed, failed, reasons = [], [], {}

    for variant in variants:
        ok, why = variant_passes_filter(variant, criteria)
        if ok:
            passed.append(variant)
        else:
            failed.append(variant)
            reasons[variant.key] = why

    return VariantFilterResult(passed=passed, failed=failed, reasons=reasons)


# ============================================================================
# Statistics
# ============================================================================

def is_transition_variant(variant: Variant) -> Optional[bool]:
    if not variant.is_snv:
        return None
    pair = {variant.reference, variant.alternative}
    return pair in ({"A", "G"}, {"C", "T"})


def calculate_variant_statistics(
    variants: Sequence[Variant],
    *,
    rare_frequency_threshold: float = 0.01,
) -> VariantStatistics:
    variants = synchronize_variant_fields(variants)
    stats = VariantStatistics()
    stats.total_variants = len(variants)
    stats.unique_variants = len({v.key for v in variants})

    qualities = []
    frequencies = []
    genes = set()
    chromosomes = set()
    consequences = Counter()

    for variant in variants:
        variant_type = variant.variant_type
        if variant_type == "SNV":
            stats.snvs += 1
            ti_tv = is_transition_variant(variant)
            if ti_tv is True:
                stats.transition_count += 1
            elif ti_tv is False:
                stats.transversion_count += 1
        elif variant_type == "MNV":
            stats.mnvs += 1
        elif variant_type == "INSERTION":
            stats.insertions += 1
            stats.indels += 1
        elif variant_type == "DELETION":
            stats.deletions += 1
            stats.indels += 1
        elif variant_type == "INDEL":
            stats.indels += 1
        elif variant_type == "COMPLEX":
            stats.complex += 1
        elif variant_type == "SYMBOLIC":
            stats.symbolic += 1
        else:
            stats.unknown += 1

        # VCF FILTER: "PASS" = passed all filters; "." = filters not
        # applied (not a failure); anything else names failed filters.
        if variant.filter_status == "PASS":
            stats.passed_filters += 1
        elif variant.filter_status in {"UNFILTERED", "UNKNOWN", "."}:
            stats.unfiltered += 1
        else:
            stats.failed_filters += 1

        if variant.quality is None:
            stats.missing_quality += 1
        else:
            qualities.append(variant.quality)

        if variant.allele_frequency is not None:
            frequencies.append(variant.allele_frequency)
            if variant.allele_frequency <= rare_frequency_threshold:
                stats.rare_variants += 1

        if variant.gene:
            genes.add(variant.gene)
        chromosomes.add(variant.chromosome)

        if variant.consequence:
            consequences[variant.consequence] += 1

    stats.genes_affected = len(genes)
    stats.chromosomes_affected = len(chromosomes)
    stats.consequence_counts = dict(consequences)

    if qualities:
        stats.mean_quality = statistics.mean(qualities)
        stats.median_quality = statistics.median(qualities)
        stats.min_quality = min(qualities)
        stats.max_quality = max(qualities)

    if frequencies:
        stats.mean_allele_frequency = statistics.mean(frequencies)
        stats.median_allele_frequency = statistics.median(frequencies)

    if stats.transversion_count:
        stats.transition_transversion_ratio = (
            stats.transition_count / stats.transversion_count
        )
    elif stats.transition_count:
        stats.transition_transversion_ratio = math.inf

    return stats


# ============================================================================
# Cohort / sample analysis
# ============================================================================

def sample_variant_summary(
    variants: Sequence[Variant],
    sample: str,
) -> SampleVariantSummary:
    summary = SampleVariantSummary(sample=sample)

    for variant in variants:
        summary.total_variant_records += 1
        gt = variant.sample_data.get(sample, {}).get("GT")
        category = genotype_category(gt)

        if category == "missing":
            summary.missing_calls += 1
        elif category == "homozygous_ref":
            summary.homozygous_ref_calls += 1
        elif category == "heterozygous":
            summary.heterozygous_calls += 1
            summary.called_variants += 1
        elif category == "homozygous_alt":
            summary.homozygous_alt_calls += 1
            summary.called_variants += 1
        elif category == "partial":
            summary.missing_calls += 1

        alleles, _ = parse_genotype(gt)
        if alleles:
            summary.allele_count += sum(
                allele is not None and allele > 0 for allele in alleles
            )

    return summary


def calculate_cohort_summary(
    variants: Sequence[Variant],
) -> CohortVariantSummary:
    samples = sorted({
        sample
        for variant in variants
        for sample in variant.sample_data
    })

    prevalence = Counter()

    for variant in variants:
        for sample, fields in variant.sample_data.items():
            gt = fields.get("GT")
            alleles, _ = parse_genotype(gt)
            if alleles and any(
                allele is not None and allele > 0
                for allele in alleles
            ):
                prevalence[variant.key] += 1

    unique_keys = {v.key for v in variants}
    shared = sum(count >= 2 for count in prevalence.values())
    private = sum(count == 1 for count in prevalence.values())

    return CohortVariantSummary(
        sample_count=len(samples),
        variant_count=len(variants),
        unique_variant_count=len(unique_keys),
        shared_variant_count=shared,
        private_variant_count=private,
        variant_prevalence=dict(prevalence),
        sample_summaries=[
            sample_variant_summary(variants, sample)
            for sample in samples
        ],
    )


# ============================================================================
# Comparison / overlap
# ============================================================================

def variant_key_set(variants: Sequence[Variant]) -> Set[Tuple[str, int, str, str]]:
    return {variant.key for variant in variants}


def compare_variant_sets(
    left: Sequence[Variant],
    right: Sequence[Variant],
) -> VariantComparison:
    left_map = {v.key: v for v in left}
    right_map = {v.key: v for v in right}

    left_keys = set(left_map)
    right_keys = set(right_map)
    shared_keys = left_keys & right_keys
    union = left_keys | right_keys

    return VariantComparison(
        only_in_left=[left_map[k] for k in sorted(left_keys - right_keys)],
        only_in_right=[right_map[k] for k in sorted(right_keys - left_keys)],
        shared=[left_map[k] for k in sorted(shared_keys)],
        left_count=len(left_keys),
        right_count=len(right_keys),
        shared_count=len(shared_keys),
        jaccard_similarity=(
            len(shared_keys) / len(union) if union else 1.0
        ),
    )


def variants_by_chromosome(
    variants: Sequence[Variant],
) -> Dict[str, List[Variant]]:
    grouped = defaultdict(list)
    for variant in variants:
        grouped[variant.chromosome].append(variant)
    return dict(grouped)


def variants_by_gene(
    variants: Sequence[Variant],
) -> Dict[str, List[Variant]]:
    grouped = defaultdict(list)
    for variant in variants:
        grouped[variant.gene or "Unknown"].append(variant)
    return dict(grouped)


def variants_by_consequence(
    variants: Sequence[Variant],
) -> Dict[str, List[Variant]]:
    grouped = defaultdict(list)
    for variant in variants:
        grouped[variant.consequence or "unknown"].append(variant)
    return dict(grouped)


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

def variants_in_interval(
    variants: Sequence[Variant],
    chromosome: str,
    start: int,
    end: int,
) -> List[Variant]:
    if start < 1 or end < start:
        raise VariantAnalysisError("Invalid interval.")
    return [
        v for v in variants
        if v.chromosome == chromosome
        and v.position <= end
        and (v.end_position or v.position) >= start
    ]


def interval_variant_density(
    variants: Sequence[Variant],
    chromosome: str,
    start: int,
    end: int,
) -> float:
    interval_length = end - start + 1
    if interval_length <= 0:
        raise VariantAnalysisError("Invalid interval.")
    count = len(variants_in_interval(variants, chromosome, start, end))
    return count / interval_length * 1000.0


# ============================================================================
# VCF serialization
# ============================================================================

def variants_to_vcf(
    variants: Sequence[Variant],
    *,
    metadata: Optional[Sequence[str]] = None,
    sample_names: Optional[Sequence[str]] = None,
) -> str:
    """Serialize variants into a valid basic VCF representation."""
    lines = [
        "##fileformat=VCFv4.3",
    ]

    if metadata:
        lines.extend(
            line if line.startswith("##") else f"##{line}"
            for line in metadata
            if line
            and not line.startswith("#CHROM")
            and line != "##fileformat=VCFv4.3"
        )

    samples = list(sample_names or sorted({
        sample
        for variant in variants
        for sample in variant.sample_data
    }))

    header = [
        "#CHROM", "POS", "ID", "REF", "ALT",
        "QUAL", "FILTER", "INFO",
    ]

    if samples:
        header.extend(["FORMAT", *samples])

    lines.append("\t".join(header))

    for variant in sorted(
        variants,
        key=lambda v: (v.chromosome, v.position, v.reference, v.alternative),
    ):
        row = [
            variant.chromosome,
            str(variant.position),
            variant.variant_id or ".",
            variant.reference,
            variant.alternative,
            "." if variant.quality is None else str(variant.quality),
            "." if variant.filter_status in (None, "", "UNFILTERED") else variant.filter_status,
            format_info_field(variant.info),
        ]

        if samples:
            format_keys = variant.format_keys or ["GT"]
            row.append(":".join(format_keys))
            for sample in samples:
                values = variant.sample_data.get(sample, {})
                row.append(
                    ":".join(str(values.get(key, ".")) for key in format_keys)
                )

        lines.append("\t".join(row))

    return "\n".join(lines) + "\n"


def write_vcf(
    variants: Sequence[Variant],
    path: str,
    *,
    metadata: Optional[Sequence[str]] = None,
    sample_names: Optional[Sequence[str]] = None,
) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(
            variants_to_vcf(
                variants,
                metadata=metadata,
                sample_names=sample_names,
            )
        )


def variants_to_tsv(variants: Sequence[Variant]) -> str:
    """Export a normalized analysis table."""
    columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "TYPE",
        "QUAL", "FILTER", "GENE", "TRANSCRIPT", "REGION",
        "CONSEQUENCE", "CODING_CHANGE", "PROTEIN_CHANGE",
        "AF", "AC", "AN",
    ]
    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)

    for v in variants:
        writer.writerow([
            v.chromosome,
            v.position,
            v.variant_id or ".",
            v.reference,
            v.alternative,
            v.variant_type,
            "." if v.quality is None else v.quality,
            v.filter_status,
            v.gene or ".",
            v.transcript or ".",
            v.region or ".",
            v.consequence or ".",
            v.coding_change or ".",
            v.protein_change or ".",
            "." if v.allele_frequency is None else v.allele_frequency,
            "." if v.allele_count is None else v.allele_count,
            "." if v.allele_number is None else v.allele_number,
        ])

    return output.getvalue()


# ============================================================================
# High-level analysis
# ============================================================================

def analyze_variants(
    variants: Sequence[Variant],
    *,
    name: str = "Variant Analysis",
    filter_criteria: Optional[VariantFilter] = None,
    annotate_info: bool = True,
    normalize: bool = False,
    reference_sequence: Optional[str] = None,
) -> VariantAnalysisResult:
    """Run the main variant-analysis workflow."""
    working = list(variants)

    for variant in working:
        issues = validate_variant(variant)
        if issues:
            variant.annotation.setdefault("validation_warnings", []).extend(issues)

    if annotate_info:
        working = annotate_from_info(working)

    if normalize:
        working = normalize_variants(
            working,
            reference_sequence=reference_sequence,
        )

    result = VariantAnalysisResult(
        name=name,
        variants=working,
        statistics=calculate_variant_statistics(working),
    )

    if filter_criteria is not None:
        result.filter_result = filter_variants(
            working,
            filter_criteria,
        )

    result.cohort_summary = calculate_cohort_summary(working)
    return result


def analyze_vcf_text(
    text: str,
    *,
    name: str = "VCF Analysis",
    filter_criteria: Optional[VariantFilter] = None,
    normalize: bool = False,
    reference_sequence: Optional[str] = None,
) -> VariantAnalysisResult:
    variants, _, _ = parse_vcf(text)
    return analyze_variants(
        variants,
        name=name,
        filter_criteria=filter_criteria,
        normalize=normalize,
        reference_sequence=reference_sequence,
    )


# ============================================================================
# Reporting
# ============================================================================

def format_variant(variant: Variant) -> str:
    quality = "." if variant.quality is None else f"{variant.quality:.2f}"
    return (
        f"{variant.chromosome}:{variant.position} "
        f"{variant.reference}>{variant.alternative} "
        f"[{variant.variant_type}] "
        f"QUAL={quality}"
    )


def build_variant_report(
    result: VariantAnalysisResult,
    *,
    max_variants: Optional[int] = 100,
) -> str:
    s = result.statistics

    lines = [
        "=" * 70,
        "SYNAPSE VARIANT ANALYSIS",
        "=" * 70,
        f"Analysis:          {result.name}",
        f"Total variants:    {s.total_variants}",
        f"Unique variants:   {s.unique_variants}",
        "",
        "VARIANT TYPES",
        "-" * 70,
        f"SNVs:              {s.snvs}",
        f"MNVs:              {s.mnvs}",
        f"Insertions:        {s.insertions}",
        f"Deletions:         {s.deletions}",
        f"Indels:            {s.indels}",
        f"Complex:           {s.complex}",
        f"Symbolic:          {s.symbolic}",
        "",
        "QUALITY",
        "-" * 70,
        f"PASS:              {s.passed_filters}",
        f"Failed filters:    {s.failed_filters}",
        f"Not filtered (.):  {s.unfiltered}",
        f"Missing QUAL:      {s.missing_quality}",
        f"Mean QUAL:         {'.' if s.mean_quality is None else f'{s.mean_quality:.2f}'}",
        f"Median QUAL:       {'.' if s.median_quality is None else f'{s.median_quality:.2f}'}",
        f"Min QUAL:          {'.' if s.min_quality is None else f'{s.min_quality:.2f}'}",
        f"Max QUAL:          {'.' if s.max_quality is None else f'{s.max_quality:.2f}'}",
        "",
        "ALLELE FREQUENCY",
        "-" * 70,
        f"Mean AF:           {'.' if s.mean_allele_frequency is None else f'{s.mean_allele_frequency:.5f}'}",
        f"Median AF:         {'.' if s.median_allele_frequency is None else f'{s.median_allele_frequency:.5f}'}",
        f"Rare variants:     {s.rare_variants}",
        "",
        "ANNOTATION",
        "-" * 70,
        f"Genes affected:    {s.genes_affected}",
        f"Chromosomes:       {s.chromosomes_affected}",
    ]

    if s.consequence_counts:
        lines.extend(["", "CONSEQUENCES", "-" * 70])
        for consequence, count in sorted(
            s.consequence_counts.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            lines.append(f"{consequence:<25} {count}")

    ratio = (
        "inf" if s.transition_transversion_ratio == math.inf
        else "N/A" if s.transition_transversion_ratio is None
        else f"{s.transition_transversion_ratio:.3f}"
    )
    lines.extend([
        "",
        "SUBSTITUTION SPECTRUM",
        "-" * 70,
        f"Transitions:       {s.transition_count}",
        f"Transversions:     {s.transversion_count}",
        f"Ti/Tv ratio:       {ratio}",
        "",
        "VARIANTS",
        "-" * 70,
    ])

    listed = result.variants if max_variants is None else result.variants[:max_variants]
    for index, variant in enumerate(listed, start=1):
        lines.append(f"{index:>4}. {format_variant(variant)}")
        if variant.gene:
            lines.append(f"       Gene: {variant.gene}")
        if variant.consequence:
            lines.append(f"       Consequence: {variant.consequence}")
        if variant.protein_change:
            lines.append(f"       Protein: {variant.protein_change}")

    if max_variants is not None and len(result.variants) > max_variants:
        lines.append(
            f"... {len(result.variants) - max_variants} additional variants omitted."
        )

    if result.filter_result is not None:
        lines.extend([
            "",
            "FILTER RESULT",
            "-" * 70,
            f"Passed:            {len(result.filter_result.passed)}",
            f"Failed:            {len(result.filter_result.failed)}",
        ])

    if result.cohort_summary:
        c = result.cohort_summary
        lines.extend([
            "",
            "COHORT",
            "-" * 70,
            f"Samples:            {c.sample_count}",
            f"Shared variants:    {c.shared_variant_count}",
            f"Private variants:   {c.private_variant_count}",
        ])

    if result.notes:
        lines.extend(["", "NOTES", "-" * 70])
        lines.extend(f"- {note}" for note in result.notes)

    lines.append("=" * 70)
    return "\n".join(lines)


def print_variant_report(
    result: VariantAnalysisResult,
    *,
    max_variants: Optional[int] = 100,
) -> None:
    print(build_variant_report(result, max_variants=max_variants))


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
    try:
        from mutation import compare_sequences
    except ImportError as exc:
        raise VariantAnalysisError(
            "mutation.py could not be imported for sequence-derived annotation."
        ) from exc

    result = compare_sequences(
        reference_sequence,
        alternative_sequence,
        reference_name="Reference",
        variant_name="Variant",
        reference_cds=reference_sequence,
        variant_cds=alternative_sequence,
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
    variant.protein_change = first.consequence
    variant.annotation["mutation_signature"] = {
        "position": first.position,
        "reference": first.reference,
        "alternative": first.alternative,
    }
    return variant


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


# ============================================================================
# Self tests
# ============================================================================

def run_variant_analysis_self_test() -> Dict[str, bool]:
    tests: Dict[str, bool] = {}

    # Basic classification.
    tests["snv_classification"] = (
        classify_variant_type("A", "G") == "SNV"
    )
    tests["mnv_classification"] = (
        classify_variant_type("AT", "GC") == "MNV"
    )
    tests["insertion_classification"] = (
        classify_variant_type("A", "ATG") == "INSERTION"
    )
    tests["deletion_classification"] = (
        classify_variant_type("ATG", "A") == "DELETION"
    )

    # INFO.
    info = parse_info_field("DP=20;AF=0.25;DB")
    tests["info_parser"] = (
        info["DP"] == 20 and
        info["AF"] == 0.25 and
        info["DB"] is True
    )

    # Genotypes.
    tests["genotype_parser"] = (
        genotype_category("0/1") == "heterozygous"
        and genotype_category("1/1") == "homozygous_alt"
        and genotype_category("0/0") == "homozygous_ref"
        and genotype_category("./.") == "missing"
    )

    # VCF parsing.
    vcf = """##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\trs1\tA\tG\t50\tPASS\tAF=0.1\tGT:DP\t0/1:20\t1/1:30
chr1\t20\t.\tC\tT,A\t40\tPASS\tDP=10\tGT\t0/0\t0/1
"""
    parsed, _, samples = parse_vcf(vcf)
    tests["vcf_parser"] = (
        len(parsed) == 3
        and samples == ["S1", "S2"]
        and parsed[0].variant_type == "SNV"
    )

    # Multi-allelic splitting.
    record = VariantRecord(
        "chr1", 100, "A", ["C", "G"], variant_id="multi"
    )
    tests["multiallelic_split"] = len(record.split()) == 2

    # Normalization.
    v = Variant("chr1", 5, "AA", "AT")
    normalized = left_normalize_variant(v)
    tests["normalization"] = normalized.normalized and normalized.position >= 1

    # Filtering.
    variants = [
        Variant("chr1", 1, "A", "G", quality=50, filter_status="PASS"),
        Variant("chr1", 2, "C", "T", quality=5, filter_status="PASS"),
    ]
    filtered = filter_variants(
        variants,
        VariantFilter(minimum_quality=20),
    )
    tests["filtering"] = (
        len(filtered.passed) == 1
        and len(filtered.failed) == 1
    )

    # Statistics.
    stats = calculate_variant_statistics(variants)
    tests["statistics"] = (
        stats.total_variants == 2
        and stats.snvs == 2
        and stats.transition_count == 2
    )

    # Set comparison.
    left = [Variant("chr1", 1, "A", "G"), Variant("chr1", 2, "C", "T")]
    right = [Variant("chr1", 2, "C", "T"), Variant("chr1", 3, "G", "A")]
    comparison = compare_variant_sets(left, right)
    tests["set_comparison"] = (
        comparison.shared_count == 1
        and comparison.jaccard_similarity == 1 / 3
    )

    # Annotation.
    annotated = annotate_variant(
        Variant("chr1", 10, "A", "G"),
        gene="GENE1",
        consequence="missense",
        allele_frequency=0.005,
    )
    tests["annotation"] = (
        annotated.gene == "GENE1"
        and annotated.consequence == "missense"
        and annotated.allele_frequency == 0.005
    )

    # Cohort.
    cohort_variants = [
        Variant(
            "chr1", 1, "A", "G",
            sample_data={
                "S1": {"GT": "0/1"},
                "S2": {"GT": "0/0"},
            },
        ),
        Variant(
            "chr1", 2, "C", "T",
            sample_data={
                "S1": {"GT": "1/1"},
                "S2": {"GT": "0/1"},
            },
        ),
    ]
    cohort = calculate_cohort_summary(cohort_variants)
    tests["cohort_analysis"] = (
        cohort.sample_count == 2
        and cohort.unique_variant_count == 2
        and cohort.shared_variant_count == 1
    )

    # Prioritization.
    prioritized = prioritize_variants([
        Variant("chr1", 1, "A", "G", consequence="synonymous"),
        Variant("chr1", 2, "A", "T", consequence="frameshift"),
    ])
    tests["prioritization"] = prioritized[0].consequence == "frameshift"

    # Interval.
    interval = variants_in_interval(
        [Variant("chr1", 10, "A", "G"), Variant("chr1", 50, "C", "T")],
        "chr1", 1, 20,
    )
    tests["interval_analysis"] = len(interval) == 1

    # Serialization.
    exported = variants_to_vcf([Variant("chr1", 1, "A", "G")])
    reparsed, _, _ = parse_vcf(exported)
    tests["vcf_roundtrip"] = (
        len(reparsed) == 1
        and reparsed[0].key == ("chr1", 1, "A", "G")
    )

    tsv = variants_to_tsv([Variant("chr1", 1, "A", "G")])
    tests["tsv_export"] = "CHROM" in tsv and "chr1" in tsv

    # Evidence.
    evidence_variant = attach_evidence(
        Variant("chr1", 1, "A", "G"),
        VariantEvidence(
            source="example",
            field="frequency",
            value=0.01,
            evidence_level="observed",
        ),
    )
    tests["evidence_container"] = (
        len(evidence_variant.annotation["evidence"]) == 1
    )

    # High-level workflow.
    workflow = analyze_variants(
        variants,
        filter_criteria=VariantFilter(minimum_quality=10),
    )
    tests["high_level_workflow"] = (
        workflow.count == 2
        and workflow.filter_result is not None
        and len(workflow.filter_result.passed) == 1
    )

    # JSON.
    tests["json_export"] = '"variants"' in workflow.to_json()

    # Validation.
    invalid = Variant("chr1", 1, "A", "G", allele_frequency=1.5)
    tests["validation"] = bool(validate_variant(invalid))

    # Integration: VCF INFO AF must feed filtering and statistics directly.
    af_vcf = """##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t100\tPASS\tAF=0.10
chr1\t200\t.\tC\tT\t100\tPASS\tAF=0.01
"""
    af_variants, _, _ = parse_vcf(af_vcf)
    af_filter = filter_variants(
        af_variants,
        VariantFilter(maximum_allele_frequency=0.05),
    )
    af_stats = calculate_variant_statistics(af_variants)
    tests["af_integration"] = (
        len(af_filter.passed) == 1
        and len(af_filter.failed) == 1
        and af_stats.mean_allele_frequency == 0.055
    )

    # Missing values are excluded when a threshold is explicitly requested.
    missing_filter = filter_variants(
        [Variant("chr1", 1, "A", "G")],
        VariantFilter(minimum_quality=10),
    )
    tests["missing_value_filtering"] = (
        len(missing_filter.passed) == 0
        and len(missing_filter.failed) == 1
    )

    # Multiallelic AF and genotype splitting must be allele-specific.
    multi_vcf = """##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t300\t.\tA\tC,G\t50\tPASS\tAF=0.1,0.2\tGT:AD\t1/2:10,5,7
"""
    multi_variants, _, _ = parse_vcf(multi_vcf)
    tests["multiallelic_allele_specific"] = (
        len(multi_variants) == 2
        and multi_variants[0].info.get("AF") == 0.1
        and multi_variants[1].info.get("AF") == 0.2
        and multi_variants[0].sample_data["S1"]["GT"] == "1/0"
        and multi_variants[1].sample_data["S1"]["GT"] == "0/1"
        and multi_variants[0].sample_data["S1"]["AD"] == "10,5"
        and multi_variants[1].sample_data["S1"]["AD"] == "10,7"
    )

    # Normalization must preserve evidence/annotation state.
    annotated_original = annotate_variant(
        Variant("chr1", 5, "AA", "AT"),
        gene="GENE1",
        consequence="missense",
    )
    attach_evidence(
        annotated_original,
        VariantEvidence(source="test", field="note", value="preserve"),
    )
    normalized_annotated = left_normalize_variant(annotated_original)
    tests["normalization_preserves_annotation"] = (
        normalized_annotated.gene == "GENE1"
        and normalized_annotated.consequence == "missense"
        and len(normalized_annotated.annotation.get("evidence", [])) == 1
    )

    # Round-trip must retain INFO AF.
    roundtrip_source = [Variant("chr1", 400, "G", "A", info={"AF": 0.25})]
    roundtrip_text = variants_to_vcf(roundtrip_source)
    roundtrip_variants, _, _ = parse_vcf(roundtrip_text)
    annotate_from_info(roundtrip_variants)
    tests["af_roundtrip"] = (
        len(roundtrip_variants) == 1
        and roundtrip_variants[0].allele_frequency == 0.25
    )

    return tests


def assert_variant_analysis_self_test() -> None:
    tests = run_variant_analysis_self_test()
    failed = [name for name, passed in tests.items() if not passed]
    if failed:
        raise AssertionError(
            "Variant analysis self-test failed: " + ", ".join(failed)
        )


# ============================================================================
# Interactive UI
# ============================================================================

def _prompt_variant() -> Variant:
    print()
    chromosome = input("Chromosome/contig: ").strip()
    position = parse_int(input("Position: ").strip())
    if position is None:
        raise VariantAnalysisError("Position must be an integer.")
    reference = input("Reference allele: ").strip()
    alternative = input("Alternative allele: ").strip()
    quality = parse_float(input("Quality (Enter for unknown): ").strip())

    variant = Variant(
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternative=alternative,
        quality=quality,
        filter_status="UNFILTERED",
    )
    return variant


def _prompt_filter() -> VariantFilter:
    print()
    print("=" * 70)
    print("VARIANT FILTER")
    print("=" * 70)

    quality_text = input("Minimum QUAL (Enter for none): ").strip()
    minimum_quality = parse_float(quality_text)

    max_af_text = input("Maximum AF (Enter for none): ").strip()
    maximum_af = parse_float(max_af_text)

    require_pass = (
        input("Require FILTER=PASS? [y/N]: ").strip().lower() == "y"
    )
    include_unknown_quality = (
        input("Include variants with unknown QUAL? [y/N]: ").strip().lower() == "y"
    )
    include_unknown_frequency = (
        input("Include variants with unknown AF? [y/N]: ").strip().lower() == "y"
    )

    return VariantFilter(
        minimum_quality=minimum_quality,
        maximum_allele_frequency=maximum_af,
        require_pass=require_pass,
        include_unknown_quality=include_unknown_quality,
        include_unknown_frequency=include_unknown_frequency,
    )


def _prompt_filename(default_name: str, extension: str) -> str:
    name = input(f"Output filename [{default_name}.{extension}]: ").strip()
    if not name:
        name = f"{default_name}.{extension}"
    elif not name.lower().endswith(f".{extension}"):
        name += f".{extension}"
    return name


def _print_filter_result(result: VariantFilterResult) -> None:
    print()
    print(f"Passed: {len(result.passed)}")
    print(f"Failed: {len(result.failed)}")
    if result.failed:
        print()
        print("FILTER FAILURES")
        print("-" * 70)
        for variant in result.failed:
            reasons = result.reasons.get(variant.key, [])
            print(f"{format_variant(variant)}")
            for reason in reasons:
                print(f"  - {reason}")


def _load_vcf_text_interactive() -> Tuple[List[Variant], List[str], List[str]]:
    print()
    path = input("VCF file path (.vcf or .vcf.gz; Enter to paste text instead): ").strip().strip('"')
    if path:
        try:
            return parse_vcf_file(path)
        except OSError as error:
            raise VariantAnalysisError(f"Could not read VCF file: {error}") from error
    print("Paste VCF text. Enter an empty line to finish.")
    lines: List[str] = []
    while True:
        line = input()
        if line == "":
            break
        lines.append(line)
    if not lines:
        raise VariantAnalysisError("No VCF text supplied.")
    return parse_vcf("\n".join(lines))


def _prompt_contig_reference(pipeline, variants: Sequence[Variant]) -> Dict[str, str]:
    """Optionally map one VCF contig to the sequence loaded in Synapse."""
    context = getattr(pipeline, "context", None)
    sequence = getattr(context, "sequence", "") if context is not None else ""
    if not sequence:
        print("No Synapse sequence loaded: variants are trimmed only (no REF check, no repeat shift).")
        return {}
    contigs = sorted({v.chromosome for v in variants})
    header = getattr(context, "header", "") or "loaded sequence"
    print()
    print(f"Loaded Synapse sequence: {header} ({len(sequence)} bp)")
    print(f"Contigs in the VCF: {', '.join(contigs)}")
    print("VCF positions must be 1-based positions on this exact sequence.")
    choice = input("Contig to map onto the loaded sequence (Enter = none): ").strip()
    if not choice:
        return {}
    if choice not in contigs:
        print(f"Contig '{choice}' is not in the VCF; no mapping used.")
        return {}
    return {choice: sequence}


def run_variant_analysis_ui(
    variants: Optional[Sequence[Variant]] = None,
    *,
    pipeline=None,
    name: str = "Variant Analysis",
) -> Optional[VariantAnalysisResult]:
    """Interactive CLI for the Variant Analysis module.

    The UI keeps the original dataset intact during filtering and normalization,
    supports real file exports, provides set comparison through VCF text input,
    and refreshes the analysis result after state-changing operations.
    """
    working = synchronize_variant_fields(list(variants or []))
    last_result: Optional[VariantAnalysisResult] = None
    last_filter_result: Optional[VariantFilterResult] = None
    metadata: List[str] = []
    sample_names: List[str] = sorted({s for v in working for s in v.sample_data})
    filtered_view: Optional[List[Variant]] = None

    def refresh_result() -> Optional[VariantAnalysisResult]:
        nonlocal last_result
        if not working:
            last_result = None
            return None
        last_result = analyze_variants(working, name=name)
        if pipeline is not None:
            _store_variant_pipeline_result(pipeline, "variant", last_result.to_dict())
        return last_result

    while True:
        print()
        print("=" * 70)
        print("VARIANT ANALYSIS")
        print("=" * 70)
        print()
        print(f"Loaded variants: {len(working)}")
        if filtered_view is not None:
            print(f"Current filter view: {len(filtered_view)} passed")
        print()
        print("INPUT")
        print("[1] Load VCF")
        print("[2] Add Variant Manually")
        print()
        print("ANALYSIS")
        print("[3] Analyze Variants")
        print("[4] Filter Variants")
        print("[5] Compare Variant Sets")
        print("[6] Variant Statistics")
        print("[7] Cohort / Sample Analysis")
        print("[8] Normalize Variants")
        print()
        print("OUTPUT")
        print("[9] Export VCF")
        print("[10] Export TSV")
        print("[11] Full Variant Report")
        print()
        print("SYSTEM")
        print("[12] Run Variant Self-Test")
        print()
        print("[0] Back")
        print("-" * 70)

        choice = input("Selection: ").strip()

        if choice == "0":
            return last_result

        if choice == "1":
            try:
                parsed, parsed_metadata, parsed_samples = _load_vcf_text_interactive()
                working = synchronize_variant_fields(parsed)
                metadata = parsed_metadata
                sample_names = parsed_samples
                filtered_view = None
                last_filter_result = None
                refresh_result()
                print(f"Loaded {len(working)} variant(s).")
            except (VariantAnalysisError, ValueError) as error:
                print(f"Error: {error}")

        elif choice == "2":
            try:
                variant = _prompt_variant()
                working.append(variant)
                working = synchronize_variant_fields(working)
                filtered_view = None
                last_filter_result = None
                refresh_result()
                print(f"Added: {format_variant(variant)}")
            except VariantAnalysisError as error:
                print(f"Error: {error}")

        elif choice == "3":
            if not working:
                print("No variants loaded.")
                continue
            try:
                result = refresh_result()
                if result is not None:
                    print()
                    print_variant_report(result)
            except VariantAnalysisError as error:
                print(f"Error: {error}")

        elif choice == "4":
            if not working:
                print("No variants loaded.")
                continue
            try:
                criteria = _prompt_filter()
                last_filter_result = filter_variants(working, criteria)
                filtered_view = list(last_filter_result.passed)
                _print_filter_result(last_filter_result)
            except VariantAnalysisError as error:
                print(f"Error: {error}")

        elif choice == "5":
            if not working:
                print("No variants loaded.")
                continue
            try:
                print()
                print("Paste the SECOND VCF set. Enter an empty line to finish.")
                lines: List[str] = []
                while True:
                    line = input()
                    if line == "":
                        break
                    lines.append(line)
                if not lines:
                    raise VariantAnalysisError("No second VCF set supplied.")
                right, _, _ = parse_vcf("\n".join(lines))
                comparison = compare_variant_sets(working, right)
                print()
                print("=" * 70)
                print("VARIANT SET COMPARISON")
                print("=" * 70)
                print(f"Left variants:     {comparison.left_count}")
                print(f"Right variants:    {comparison.right_count}")
                print(f"Shared variants:   {comparison.shared_count}")
                print(f"Jaccard similarity: {comparison.jaccard_similarity:.4f}")
                print()
                print("Only in current set:")
                for v in comparison.only_in_left:
                    print(f"  - {format_variant(v)}")
                print("Only in second set:")
                for v in comparison.only_in_right:
                    print(f"  - {format_variant(v)}")
            except (VariantAnalysisError, ValueError) as error:
                print(f"Error: {error}")

        elif choice == "6":
            if not working:
                print("No variants loaded.")
                continue
            stats = calculate_variant_statistics(working)
            print()
            print("=" * 70)
            print("VARIANT STATISTICS")
            print("=" * 70)
            print(f"Total:              {stats.total_variants}")
            print(f"Unique:             {stats.unique_variants}")
            print(f"SNVs:               {stats.snvs}")
            print(f"MNVs:               {stats.mnvs}")
            print(f"Insertions:         {stats.insertions}")
            print(f"Deletions:          {stats.deletions}")
            print(f"Indels:             {stats.indels}")
            print(f"Missing QUAL:       {stats.missing_quality}")
            print(f"Mean QUAL:          {'.' if stats.mean_quality is None else f'{stats.mean_quality:.2f}'}")
            print(f"Mean AF:            {'.' if stats.mean_allele_frequency is None else f'{stats.mean_allele_frequency:.5f}'}")
            print(f"Median AF:          {'.' if stats.median_allele_frequency is None else f'{stats.median_allele_frequency:.5f}'}")
            print(f"Rare variants:      {stats.rare_variants}")
            print(f"Genes affected:     {stats.genes_affected}")
            print(f"Chromosomes:        {stats.chromosomes_affected}")
            print(f"Transitions:        {stats.transition_count}")
            print(f"Transversions:      {stats.transversion_count}")
            ratio = stats.transition_transversion_ratio
            print(f"Ti/Tv ratio:        {'inf' if ratio == math.inf else '.' if ratio is None else f'{ratio:.3f}'}")

        elif choice == "7":
            if not working:
                print("No variants loaded.")
                continue
            cohort = calculate_cohort_summary(working)
            print()
            print("=" * 70)
            print("COHORT / SAMPLE ANALYSIS")
            print("=" * 70)
            print(f"Samples:            {cohort.sample_count}")
            print(f"Shared variants:    {cohort.shared_variant_count}")
            print(f"Private variants:   {cohort.private_variant_count}")
            for sample in cohort.sample_summaries:
                print(
                    f"{sample.sample}: called={sample.called_variants}, "
                    f"het={sample.heterozygous_calls}, "
                    f"hom_alt={sample.homozygous_alt_calls}, "
                    f"missing={sample.missing_calls}, "
                    f"alt_alleles={sample.allele_count}"
                )

        elif choice == "8":
            if not working:
                print("No variants loaded.")
                continue
            try:
                contig_map = _prompt_contig_reference(pipeline, working)
                if contig_map:
                    report = check_reference_alleles(working, contig_map)
                    print(
                        f"REF check against the loaded sequence: {report['checked']} checked, "
                        f"{report['mismatched']} mismatch(es), {report['not_checked']} not checked."
                    )
                working = normalize_variants(working, contig_references=contig_map or {})
                working = synchronize_variant_fields(working)
                filtered_view = None
                last_filter_result = None
                refresh_result()
                print(f"Normalized {len(working)} variant(s).")
            except VariantAnalysisError as error:
                print(f"Error: {error}")

        elif choice == "9":
            if not working:
                print("No variants loaded.")
                continue
            try:
                filename = _prompt_filename("synapse_variants", "vcf")
                write_vcf(working, filename, metadata=metadata, sample_names=sample_names)
                print(f"VCF exported: {filename}")
            except OSError as error:
                print(f"Error writing VCF: {error}")

        elif choice == "10":
            if not working:
                print("No variants loaded.")
                continue
            try:
                filename = _prompt_filename("synapse_variants", "tsv")
                with open(filename, "w", encoding="utf-8", newline="") as handle:
                    handle.write(variants_to_tsv(working))
                print(f"TSV exported: {filename}")
            except OSError as error:
                print(f"Error writing TSV: {error}")

        elif choice == "11":
            if not working:
                print("No variants loaded.")
                continue
            result = refresh_result()
            if result is not None:
                print()
                print_variant_report(result)

        elif choice == "12":
            tests = run_variant_analysis_self_test()
            print()
            print("=" * 70)
            print("VARIANT ANALYSIS SELF-TEST")
            print("=" * 70)
            for test_name, passed in tests.items():
                print(f"[{'PASS' if passed else 'FAIL'}] {test_name}")
            print("-" * 70)
            print(f"{sum(tests.values())}/{len(tests)} tests passed.")

        else:
            print("Invalid selection.")


def _store_variant_pipeline_result(pipeline, key: str, value: object) -> None:
    """Store a result without making the pipeline a hard dependency."""
    import logging
    logger = logging.getLogger("synapse.variant")
    if pipeline is None:
        return
    try:
        if hasattr(pipeline, "store") and callable(pipeline.store):
            pipeline.store(key, value)
        elif hasattr(pipeline, "results") and isinstance(pipeline.results, dict):
            pipeline.results[key] = value
        elif isinstance(pipeline, dict):
            pipeline[key] = value
        else:
            logger.warning("Pipeline object has no result store; %s result not stored.", key)
    except Exception:
        # Analysis success must not be lost because pipeline integration
        # differs, but the failure is logged (audit Phase V).
        logger.exception("Could not store %s result in the pipeline.", key)


# ============================================================================
# Public API
# ============================================================================

__all__ = [
    "VariantAnalysisError",
    "Variant",
    "VariantRecord",
    "VariantFilter",
    "VariantFilterResult",
    "VariantStatistics",
    "SampleVariantSummary",
    "CohortVariantSummary",
    "VariantComparison",
    "VariantAnalysisResult",
    "VariantEvidence",
    "parse_float",
    "parse_int",
    "validate_allele",
    "validate_variant",
    "classify_variant_type",
    "parse_info_field",
    "format_info_field",
    "parse_genotype",
    "genotype_category",
    "parse_vcf",
    "parse_vcf_file",
    "trim_common_bases",
    "left_normalize_variant",
    "normalize_variants",
    "annotate_variant",
    "annotate_from_info",
    "synchronize_variant_fields",
    "variant_passes_filter",
    "filter_variants",
    "is_transition_variant",
    "calculate_variant_statistics",
    "sample_variant_summary",
    "calculate_cohort_summary",
    "variant_key_set",
    "compare_variant_sets",
    "variants_by_chromosome",
    "variants_by_gene",
    "variants_by_consequence",
    "consequence_priority",
    "prioritize_variants",
    "variants_in_interval",
    "interval_variant_density",
    "variants_to_vcf",
    "write_vcf",
    "variants_to_tsv",
    "analyze_variants",
    "analyze_vcf_text",
    "format_variant",
    "build_variant_report",
    "print_variant_report",
    "annotate_from_mutation_analysis",
    "attach_evidence",
    "highest_priority_consequence",
    "run_variant_analysis_self_test",
    "assert_variant_analysis_self_test",
    "run_variant_analysis_ui",
]


if __name__ == "__main__":
    tests = run_variant_analysis_self_test()
    passed = sum(tests.values())
    print("=" * 70)
    print("SYNAPSE VARIANT ANALYSIS SELF-TEST")
    print("=" * 70)
    for name, success in tests.items():
        print(f"[{'PASS' if success else 'FAIL'}] {name}")
    print("-" * 70)
    print(f"{passed}/{len(tests)} tests passed.")
    if passed != len(tests):
        raise SystemExit(1)
