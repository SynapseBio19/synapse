"""Load sequence files with provenance (format, path, SHA-256, topology).

Supported formats:

* FASTA (``.fa``, ``.fasta``, ``.fna``, ``.ffn``, ``.txt``): parsed by the
  existing ``fileio.fasta.load_fasta``. FASTA has no topology field; a header
  token ``topology=circular`` or ``[topology=circular]`` (NCBI-style
  modifier) marks a circular molecule, otherwise records are linear.
* GenBank (``.gb``, ``.gbk``, ``.genbank``): ``fileio.genbank.load_genbank``;
  topology comes from the LOCUS line.

The format is detected from the first non-empty line, not only the
extension.
"""

from __future__ import annotations

import os
import re
from typing import Any, Dict, List

from core.models import sha256_file
from fileio.fasta import load_fasta
from fileio.genbank import load_genbank

_TOPOLOGY_TOKEN = re.compile(r"\[?topology\s*=\s*(circular|linear)\]?", re.IGNORECASE)


def detect_format(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                return "fasta"
            if stripped.startswith("LOCUS"):
                return "genbank"
            break
    raise ValueError(
        "Unrecognised sequence file: expected FASTA ('>' header) or GenBank ('LOCUS')."
    )


def fasta_topology(header: str) -> str:
    match = _TOPOLOGY_TOKEN.search(header or "")
    return match.group(1).lower() if match else "linear"


def load_sequences(path: str) -> List[Dict[str, Any]]:
    """Load all records of a FASTA or GenBank file, with provenance fields."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Sequence file not found: {path}")

    file_format = detect_format(path)
    if file_format == "fasta":
        records = load_fasta(path)
        for record in records:
            record["topology"] = fasta_topology(record["header"])
    else:
        records = load_genbank(path)

    digest = sha256_file(path)
    absolute = os.path.abspath(path)
    for index, record in enumerate(records, start=1):
        record["source_path"] = absolute
        record["source_format"] = file_format
        record["source_sha256"] = digest
        record["record_index"] = index
    return records
