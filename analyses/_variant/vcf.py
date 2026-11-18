"""VCF parsing: header numbers, INFO and FORMAT fields, genotypes, multi-allelic split."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from analyses._variant.models import (
    DNA_BASES,
    GT_ALLELE_SEPARATOR_RE,
    VALID_FILTER_STATUS,
    Variant,
    VariantAnalysisError,
)


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
    genotype: Optional[str],
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
    skipped: Optional[List[str]] = None,
) -> Tuple[List[Variant], List[str], List[str]]:
    """
    Parse VCF text.

    Returns:
        variants, metadata_lines, sample_names

    A record without an ALT allele (ALT ".", a reference-only site such as a
    gVCF reference block) is not a variant and yields no Variant when
    records are split. Pass a list as ``skipped`` to receive one
    "CHROM:POS REF>." entry per such record, so that callers can report
    them instead of losing them silently (Master Prompt 10.8).
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
            split_variants = record.split()
            if not split_variants and skipped is not None:
                skipped.append(f"{chrom}:{position} {ref}>{alt_text}")
            variants.extend(split_variants)
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


def parse_vcf_file(path: str, *, split_multiallelic: bool = True,
                   skipped: Optional[List[str]] = None):
    """Parse a VCF file; ``.gz``/BGZF-compressed files are read transparently.

    ``skipped`` collects the reference-only records, as in parse_vcf().
    """
    import gzip

    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return parse_vcf(handle.read(), split_multiallelic=split_multiallelic, skipped=skipped)
    with open(path, "r", encoding="utf-8") as handle:
        return parse_vcf(handle.read(), split_multiallelic=split_multiallelic, skipped=skipped)
