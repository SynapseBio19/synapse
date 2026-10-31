"""Minimal GenBank flat-file reader (standard library only).

Reads, per record: LOCUS name, declared length, molecule type and
topology, DEFINITION, ACCESSION, the ORIGIN sequence and the FEATURES table
(feature key, raw location string, qualifiers, and intervals for simple
locations).

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
    """Return one dict per record: header, sequence, topology, molecule type, features, notes."""
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
        molecule_type = None
        for i, token in enumerate(tokens):
            if token in {"bp", "aa"} and i > 0 and tokens[i - 1].isdigit():
                declared_length = int(tokens[i - 1])
                # GenPept records give the length in "aa" and no molecule
                # token; nucleotide records name it after "bp" (DNA, mRNA, ss-RNA, ...).
                if token == "aa":
                    molecule_type = "protein"
                elif i + 1 < len(tokens) and ("DNA" in tokens[i + 1].upper() or "RNA" in tokens[i + 1].upper()):
                    molecule_type = tokens[i + 1]
        topology = "circular" if "circular" in tokens else "linear"

        definition_parts: List[str] = []
        accession = None
        feature_lines: List[str] = []
        sequence_parts: List[str] = []
        dropped_origin_symbols: set = set()
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
                # Digits are ORIGIN line numbers and whitespace separates the
                # blocks, so both are dropped by design. Any other symbol would
                # be sequence content this reader cannot keep (audit S21).
                dropped_origin_symbols.update(
                    ch for ch in line if not ch.isalnum() and not ch.isspace()
                )

        sequence = "".join(sequence_parts).upper()
        notes: List[str] = []
        if not sequence:
            raise GenBankFormatError(f"Record '{name}' has no ORIGIN sequence.")
        if declared_length is not None and declared_length != len(sequence):
            notes.append(
                f"LOCUS declares {declared_length} bp but ORIGIN contains {len(sequence)} bp."
            )
        if dropped_origin_symbols:
            shown = ", ".join(repr(symbol) for symbol in sorted(dropped_origin_symbols))
            notes.append(
                f"ORIGIN contained symbols that are not sequence letters and were not "
                f"read: {shown}. Only letters are kept, so the stored sequence is "
                f"shorter than the file suggests."
            )

        features = _parse_features(feature_lines)
        # Locations in forms this reader does not parse (order(), remote
        # accession references, between-base sites) keep only their raw string
        # and have no interval. Say so per record instead of silently losing
        # the feature's position (audit S10).
        unparsed: Dict[str, int] = {}
        for feature in features:
            if not feature["location_parsed"]["parsed"]:
                unparsed[feature["key"]] = unparsed.get(feature["key"], 0) + 1
        if unparsed:
            shown = ", ".join(
                key if count == 1 else f"{key} x{count}"
                for key, count in sorted(unparsed.items())
            )
            notes.append(
                f"{sum(unparsed.values())} feature location(s) use a form this reader "
                f"does not parse; they keep their raw location string and have no "
                f"interval: {shown}."
            )

        definition = " ".join(definition_parts).strip()
        header = f"{accession or name} {definition}".strip()
        records.append({
            "header": header,
            "sequence": sequence,
            "topology": topology,
            "molecule_type": molecule_type,
            "locus_name": name,
            "accession": accession,
            "features": features,
            "notes": notes,
        })

    if not records:
        raise GenBankFormatError("No GenBank records found.")
    return records
