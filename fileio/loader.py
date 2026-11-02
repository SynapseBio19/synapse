"""Load sequence files with provenance (format, path, SHA-256, topology).

Supported formats:

* FASTA (``.fa``, ``.fasta``, ``.fna``, ``.ffn``, ``.txt``): parsed by the
  existing ``fileio.fasta.load_fasta``. FASTA has no topology field; a header
  token ``topology=circular`` or ``[topology=circular]`` (NCBI-style
  modifier) marks a circular molecule, otherwise records are linear.
* GenBank (``.gb``, ``.gbk``, ``.genbank``): ``fileio.genbank.load_genbank``;
  topology comes from the LOCUS line.

Sequence type: every record is classified with
``core.alphabets.classify_sequence`` (``alphabet``,
``sequence_classification``). A declared molecule type makes it CONFIRMED:
the GenBank LOCUS molecule (``aa`` for protein records, ``DNA``/``mRNA``/...
for nucleotide records) or the NCBI FASTA modifier ``[moltype=...]`` (e.g.
``[moltype=mRNA]``, ``[moltype=genomic DNA]``); ``[moltype=protein]`` is a
Synapse extension. NCBI assumes genomic DNA when the modifier is absent;
Synapse does not, and classifies the letters instead.

The format is detected from the first non-empty line, not only the
extension.
"""

from __future__ import annotations

import os
import re
from typing import Any, Dict, List, Optional

from core.alphabets import INVALID, classify_sequence, declared_type_from_molecule
from core.models import sha256_file
from fileio.fasta import load_fasta
from fileio.genbank import load_genbank

_TOPOLOGY_TOKEN = re.compile(r"\[?topology\s*=\s*(circular|linear)\]?", re.IGNORECASE)
_MOLTYPE_TOKEN = re.compile(r"\[moltype=([^\]]+)\]", re.IGNORECASE)


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


def fasta_molecule_type(header: str) -> Optional[str]:
    """Value of an NCBI ``[moltype=...]`` header modifier, or None."""
    match = _MOLTYPE_TOKEN.search(header or "")
    return match.group(1).strip() if match else None


def classify_record(record: Dict[str, Any]) -> None:
    """Set ``declared_type``, ``alphabet`` and ``sequence_classification`` on *record*.

    A declared molecule type the letters contradict is reported in
    ``notes``; the record is kept so the user sees why it is rejected.
    """
    molecule = record.get("molecule_type")
    declared = declared_type_from_molecule(molecule)
    result = classify_sequence(record.get("sequence", ""), declared)
    record["declared_type"] = declared
    record["alphabet"] = result.sequence_type
    record["sequence_classification"] = result.as_dict()
    if declared is not None and result.evidence == INVALID:
        record.setdefault("notes", []).append(
            f"Declared molecule type {molecule!r} does not fit the sequence: {result.reasons[0]}."
        )


def load_sequences(path: str) -> List[Dict[str, Any]]:
    """Load all records of a FASTA or GenBank file, with provenance fields."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Sequence file not found: {path}")

    file_format = detect_format(path)
    if file_format == "fasta":
        records = load_fasta(path)
        for record in records:
            record["topology"] = fasta_topology(record["header"])
            record["molecule_type"] = fasta_molecule_type(record["header"])
    else:
        records = load_genbank(path)

    digest = sha256_file(path)
    absolute = os.path.abspath(path)
    for index, record in enumerate(records, start=1):
        record["source_path"] = absolute
        record["source_format"] = file_format
        record["source_sha256"] = digest
        record["record_index"] = index
        classify_record(record)
    return records
