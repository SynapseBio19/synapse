"""Shared helpers for the Synapse test suite.

Importing this module:

* puts the project root on ``sys.path``;
* blocks all outgoing network connections, so no test can silently depend
  on NCBI, Europe PMC, OpenAlex or Crossref being reachable;
* replaces ``input()`` with a function that fails loudly, so a test that
  accidentally enters an interactive menu fails instead of hanging.
"""

from __future__ import annotations

import builtins
import json
import math
import os
import socket
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TESTS_DIR = Path(__file__).resolve().parent
FIXTURES_DIR = TESTS_DIR / "fixtures"
GOLDEN_DIR = TESTS_DIR / "golden"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


class NetworkBlockedError(OSError):
    """Raised when a test tries to open a network connection."""


def _blocked(*_args, **_kwargs):
    raise NetworkBlockedError("Network access is blocked in the Synapse test suite.")


socket.socket.connect = _blocked  # type: ignore[assignment]
socket.create_connection = _blocked  # type: ignore[assignment]
socket.getaddrinfo = _blocked  # type: ignore[assignment]


def _no_input(*_args, **_kwargs):
    raise RuntimeError("input() was called during an automated test.")


builtins.input = _no_input  # type: ignore[assignment]


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    """Minimal FASTA reader used only by tests (independent of fileio)."""
    records: list[tuple[str, str]] = []
    header = None
    parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header, parts = line[1:].strip(), []
        else:
            parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def reverse_complement(sequence: str) -> str:
    """Independent reverse complement for A/C/G/T (test oracle)."""
    table = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(table[b] for b in reversed(sequence))


def deterministic_dna(length: int, seed: int = 12345) -> str:
    """Reproducible pseudo-random DNA (LCG; identical on every Python version)."""
    state = seed & 0xFFFFFFFF
    out = []
    for _ in range(length):
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        out.append("ACGT"[(state >> 16) & 3])
    return "".join(out)


def normalize_for_json(value):
    """Convert analysis output into stable, JSON-comparable data."""
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return round(value, 6)
    if isinstance(value, dict):
        return {str(k): normalize_for_json(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))}
    if isinstance(value, (list, tuple)):
        return [normalize_for_json(v) for v in value]
    if isinstance(value, set):
        return sorted(normalize_for_json(v) for v in value)
    return value


def golden_path(name: str) -> Path:
    return GOLDEN_DIR / f"{name}.json"


def update_golden_requested() -> bool:
    return os.environ.get("SYNAPSE_UPDATE_GOLDEN", "") == "1"


def load_golden(name: str):
    return json.loads(golden_path(name).read_text(encoding="utf-8"))


def write_golden(name: str, data) -> None:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    golden_path(name).write_text(
        json.dumps(normalize_for_json(data), indent=1, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
