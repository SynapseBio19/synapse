"""Shared, coordinate-safe adapter for handing a selected ORF to other modules.

Background (audit P0-2, 2026-10-29): CRISPR and BLAST each had their own ORF
adapter. Both read the Synapse ORF field ``start`` (0-based) as if it were
1-based, and both ignored the strand. For a minus-strand ORF the ORF sequence
is the reverse complement of the genomic region, so every CRISPR guide
coordinate was mirrored.

Synapse ORF dictionaries (``analyses.orf.find_orfs``) use:

* ``start`` / ``stop``: 0-based, inclusive, forward-strand coordinates with
  ``start <= stop`` for both strands;
* ``strand``: "+" or "-";
* ``sequence``: the ORF in its own 5'->3' coding orientation (for "-" this is
  the reverse complement of ``dna[start:stop + 1]``).

Generic ORF-like objects from elsewhere (tests, external callers) may carry a
1-based ``start``; they are recognized by the absence of Synapse fields and
can state ``coordinate_base`` explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Optional, Sequence

_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}


def reverse_complement(sequence: str) -> str:
    return "".join(_COMPLEMENT.get(base, "N") for base in reversed(sequence.upper()))


@dataclass(frozen=True)
class OrfRegion:
    """An ORF located on a parent DNA sequence.

    ``start_1based`` / ``end_1based`` are inclusive forward-strand positions on
    the parent sequence (None when the parent position is unknown).
    """

    name: str
    coding_sequence: str
    strand: str
    start_1based: Optional[int]
    end_1based: Optional[int]
    forward_region: Optional[str]
    validated: bool
    note: str = ""

    @property
    def length(self) -> int:
        return len(self.coding_sequence)


def _get(orf: Any, keys: Sequence[str], default: Any = None) -> Any:
    if isinstance(orf, Mapping):
        for key in keys:
            if key in orf and orf[key] is not None:
                return orf[key]
    for key in keys:
        if hasattr(orf, key) and getattr(orf, key) is not None:
            return getattr(orf, key)
    return default


def is_synapse_orf(orf: Any) -> bool:
    """True for ORF dictionaries produced by ``analyses.orf.find_orfs``."""
    return all(_get(orf, (key,)) is not None for key in ("strand", "stop", "frame"))


def orf_display_name(orf: Any) -> str:
    explicit = _get(orf, ("name", "id", "orf_id", "label"))
    if explicit:
        return str(explicit)
    rank = _get(orf, ("rank",))
    if rank is not None:
        try:
            return f"ORF-{int(rank):03d}"
        except (TypeError, ValueError):
            pass
    return "Selected ORF"


def resolve_orf_region(orf: Any, parent_sequence: Optional[str] = None) -> Optional[OrfRegion]:
    """Locate an ORF on its parent sequence with explicit strand handling.

    Returns None when the ORF carries no sequence.
    """
    if orf is None:
        return None

    coding = "".join(str(_get(orf, ("sequence", "dna", "orf_sequence", "seq"), "")).split()).upper()
    if not coding:
        return None

    name = orf_display_name(orf)
    strand = str(_get(orf, ("strand",), "+")).strip() or "+"
    if strand not in {"+", "-"}:
        strand = "+"

    base = _get(orf, ("coordinate_base",))
    raw_start = _get(orf, ("start", "start_position", "genomic_start", "position"))

    if raw_start is None:
        return OrfRegion(name, coding, strand, None, None, None, False,
                         "ORF has no parent coordinates; positions are relative to the ORF.")

    try:
        raw_start = int(raw_start)
    except (TypeError, ValueError):
        return OrfRegion(name, coding, strand, None, None, None, False,
                         "ORF start is not an integer; positions are relative to the ORF.")

    if base is None:
        base = 0 if is_synapse_orf(orf) else 1
    start_1based = raw_start + 1 if int(base) == 0 else raw_start
    end_1based = start_1based + len(coding) - 1

    forward_region = None
    validated = False
    note = ""

    if parent_sequence:
        parent = "".join(str(parent_sequence).split()).upper()
        if 1 <= start_1based and end_1based <= len(parent):
            forward_region = parent[start_1based - 1:end_1based]
            expected = coding if strand == "+" else reverse_complement(coding)
            validated = forward_region == expected
            if not validated:
                note = ("ORF sequence does not match the parent sequence at the stated "
                        "coordinates; parent coordinates are not used.")
                forward_region = None
        else:
            note = "ORF coordinates lie outside the parent sequence."

    if forward_region is None and strand == "+":
        forward_region = coding if not parent_sequence else None

    return OrfRegion(name, coding, strand, start_1based, end_1based,
                     forward_region, validated, note)
