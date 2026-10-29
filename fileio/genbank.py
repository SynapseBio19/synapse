"""Minimal GenBank flat-file reader (standard library only).

Reads, per record: LOCUS name, declared length and topology, DEFINITION,
ACCESSION, the ORIGIN sequence and the FEATURES table (feature key, raw
location string, qualifiers, and intervals for simple locations).

Scope and limitations:
* Locations are parsed for the forms ``a..b``, ``complement(a..b)``,
  ``join(...)`` and ``complement(join(...))`` with plain integers. Fuzzy
  positions (``<1..>200``) are accepted with the markers stripped and the
  feature flagged ``partial``. Remote references (``AB123.1:10..20``),
  ``order()`` and between-base sites (``10^11``) keep only the raw string.
* No validation against the INSDC feature table specification beyond this.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

_LOCATION_SPAN = re.compile(r"^<?(\d+)\.\.>?(\d+)$")
_LOCATION_SINGLE = re.compile(r"^<?(\d+)>?$")


class GenBankFormatError(ValueError):
    """Raised for files that are not readable GenBank records."""


def _parse_spans(text: str) -> Optional[List[Tuple[int, int]]]:
    spans = []
    for part in text.split(","):
        part = part.strip()
        match = _LOCATION_SPAN.match(part) or _LOCATION_SINGLE.match(part)
        if not match:
            return None
        start = int(match.group(1))
        end = int(match.group(2)) if match.lastindex and match.lastindex >= 2 else start
        spans.append((start, end))
    return spans


def parse_location(location: str) -> Dict[str, Any]:
    """Parse a GenBank location into 1-based inclusive spans + strand."""
    raw = location.replace(" ", "")
    strand = "+"
    inner = raw
    if inner.startswith("complement(") and inner.endswith(")"):
        strand = "-"
        inner = inner[len("complement("):-1]
    if inner.startswith("join(") and inner.endswith(")"):
        inner = inner[len("join("):-1]
    spans = _parse_spans(inner) if "complement(" not in inner else None
    return {
        "raw": location,
        "strand": strand,
        "spans_1based": spans,
        "partial": "<" in raw or ">" in raw,
        "parsed": spans is not None,
    }


def _parse_features(lines: List[str]) -> List[Dict[str, Any]]:
    features: List[Dict[str, Any]] = []
    current: Optional[Dict[str, Any]] = None
    qualifier: Optional[str] = None
    for line in lines:
        if len(line) < 21:
            continue
        key = line[5:21].strip()
        body = line[21:].rstrip()
        if key:
            current = {"key": key, "location": body, "qualifiers": {}}
            features.append(current)
            qualifier = None
        elif current is not None:
            if body.startswith("/"):
                name, _, value = body[1:].partition("=")
                current["qualifiers"][name] = value.strip('"')
                qualifier = name
            elif qualifier is None:
                current["location"] += body.strip()
            else:
                previous = current["qualifiers"][qualifier]
                joiner = "" if qualifier == "translation" else " "
                current["qualifiers"][qualifier] = (previous + joiner + body.strip().strip('"')).strip()
    for feature in features:
        feature["location_parsed"] = parse_location(feature["location"])
    return features


def load_genbank(path: str) -> List[Dict[str, Any]]:
    """Return one dict per record: header, sequence, topology, features, notes."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        text = handle.read()

    if not text.lstrip().startswith("LOCUS"):
        raise GenBankFormatError("GenBank file must start with a LOCUS line.")

    records: List[Dict[str, Any]] = []
    for block in re.split(r"^//\s*$", text, flags=re.MULTILINE):
        if not block.strip():
            continue
        lines = block.strip("\n").splitlines()
        locus_line = next((l for l in lines if l.startswith("LOCUS")), None)
        if locus_line is None:
            raise GenBankFormatError("Record without LOCUS line.")

        tokens = locus_line.split()
        name = tokens[1] if len(tokens) > 1 else "unnamed"
        declared_length = None
        for i, token in enumerate(tokens):
            if token in {"bp", "aa"} and i > 0 and tokens[i - 1].isdigit():
                declared_length = int(tokens[i - 1])
        topology = "circular" if "circular" in tokens else "linear"

        definition_parts: List[str] = []
        accession = None
        feature_lines: List[str] = []
        sequence_parts: List[str] = []
        section = None
        for line in lines:
            if line.startswith("DEFINITION"):
                section = "definition"
                definition_parts.append(line[12:].strip())
            elif line.startswith("ACCESSION"):
                section = None
                accession = line[12:].split()[0] if line[12:].split() else None
            elif line.startswith("FEATURES"):
                section = "features"
            elif line.startswith("ORIGIN"):
                section = "origin"
            elif line[:1].isalpha():
                section = None
            elif section == "definition" and line.startswith(" "):
                definition_parts.append(line.strip())
            elif section == "features":
                feature_lines.append(line)
            elif section == "origin":
                sequence_parts.append("".join(ch for ch in line if ch.isalpha()))

        sequence = "".join(sequence_parts).upper()
        notes: List[str] = []
        if not sequence:
            raise GenBankFormatError(f"Record '{name}' has no ORIGIN sequence.")
        if declared_length is not None and declared_length != len(sequence):
            notes.append(
                f"LOCUS declares {declared_length} bp but ORIGIN contains {len(sequence)} bp."
            )
        definition = " ".join(definition_parts).strip()
        header = f"{accession or name} {definition}".strip()
        records.append({
            "header": header,
            "sequence": sequence,
            "topology": topology,
            "locus_name": name,
            "accession": accession,
            "features": _parse_features(feature_lines),
            "notes": notes,
        })

    if not records:
        raise GenBankFormatError("No GenBank records found.")
    return records
