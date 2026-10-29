"""Shared nucleotide utilities and NCBI genetic codes.

One source of truth for reverse complements, IUPAC symbols and translation
tables. Existing analysis modules still carry their own copies; they are
cross-checked against this module in tests/test_seqtools.py and will be
migrated without behavior change (see docs/CHANGES.md).

Genetic codes are written in the compact NCBI format (as in NCBI's
``gc.prt``): 64 amino-acid letters and 64 start flags in TCAG x TCAG x TCAG
codon order. Only the tables used by Synapse are included.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Dict, FrozenSet, Iterable, List, Optional

DNA_BASES = frozenset("ACGT")

IUPAC_DNA: Dict[str, FrozenSet[str]] = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"), "T": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("CG"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
}

IUPAC_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}


def normalize_dna(sequence: str) -> str:
    """Remove whitespace, uppercase, convert U to T. Does not drop other symbols."""
    return "".join(str(sequence or "").split()).upper().replace("U", "T")


def invalid_symbols(sequence: str) -> List[int]:
    """1-based positions of symbols that are not IUPAC DNA."""
    return [i for i, base in enumerate(sequence, start=1) if base not in IUPAC_DNA]


def complement(sequence: str) -> str:
    try:
        return "".join(IUPAC_COMPLEMENT[base] for base in sequence)
    except KeyError as exc:
        raise ValueError(f"Cannot complement non-IUPAC symbol {exc.args[0]!r}.") from None


def reverse_complement(sequence: str) -> str:
    return complement(sequence)[::-1]


def gc_fraction(sequence: str) -> Optional[float]:
    """GC fraction over unambiguous A/C/G/T bases; None if there are none."""
    counted = [b for b in sequence if b in DNA_BASES]
    if not counted:
        return None
    return sum(b in "GC" for b in counted) / len(counted)


# ---------------------------------------------------------------------------
# NCBI genetic codes
# ---------------------------------------------------------------------------

_BASE1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
_BASE2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
_BASE3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

_NCBI_TABLES = {
    1: ("Standard",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M------**--*----M---------------M----------------------------"),
    2: ("Vertebrate Mitochondrial",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "----------**--------------------MMMM----------**---M------------"),
    4: ("Mold, Protozoan, Coelenterate Mitochondrial; Mycoplasma/Spiroplasma",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM------**-------M------------MMMM---------------M------------"),
    11: ("Bacterial, Archaeal and Plant Plastid",
         "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M------**--*----M------------MMMM---------------M------------"),
}


@dataclass(frozen=True)
class GeneticCode:
    table_id: int
    name: str
    codons: Dict[str, str]
    start_codons: FrozenSet[str]
    stop_codons: FrozenSet[str]

    def translate_codon(self, codon: str) -> str:
        """Amino acid for a codon; '*' for stop, 'X' if it cannot be resolved.

        Codons with IUPAC symbols resolve when every expansion gives the same
        amino acid (e.g. GCN -> A); otherwise 'X'.
        """
        codon = codon.upper().replace("U", "T")
        if codon in self.codons:
            return self.codons[codon]
        if len(codon) != 3 or any(b not in IUPAC_DNA for b in codon):
            return "X"
        meanings = {self.codons["".join(c)] for c in product(*(sorted(IUPAC_DNA[b]) for b in codon))}
        return meanings.pop() if len(meanings) == 1 else "X"

    def translate(self, sequence: str, *, initiator_start: bool = False,
                  to_stop: bool = True) -> str:
        """Translate frame +1. initiator_start decodes a start codon at position 1 as M."""
        sequence = normalize_dna(sequence)
        protein: List[str] = []
        for index in range(0, len(sequence) - 2, 3):
            codon = sequence[index:index + 3]
            if index == 0 and initiator_start and codon in self.start_codons:
                protein.append("M")
                continue
            amino_acid = self.translate_codon(codon)
            if amino_acid == "*" and to_stop:
                break
            protein.append(amino_acid)
        return "".join(protein)


def _build(table_id: int) -> GeneticCode:
    name, amino_acids, starts = _NCBI_TABLES[table_id]
    codons = {b1 + b2 + b3: aa for b1, b2, b3, aa in zip(_BASE1, _BASE2, _BASE3, amino_acids)}
    start_codons = frozenset(c for c, flag in zip(codons, starts) if flag == "M")
    stop_codons = frozenset(c for c, aa in codons.items() if aa == "*")
    return GeneticCode(table_id, name, codons, start_codons, stop_codons)


GENETIC_CODES: Dict[int, GeneticCode] = {table_id: _build(table_id) for table_id in _NCBI_TABLES}


def get_genetic_code(table_id: int = 1) -> GeneticCode:
    try:
        return GENETIC_CODES[int(table_id)]
    except (KeyError, ValueError, TypeError):
        raise ValueError(
            f"Genetic code {table_id!r} is not available. "
            f"Supported NCBI tables: {sorted(GENETIC_CODES)}."
        ) from None


def available_genetic_codes() -> Iterable[GeneticCode]:
    return [GENETIC_CODES[k] for k in sorted(GENETIC_CODES)]


def codon_dict(table_id: int = 1, *, rna: bool = False, stop_symbol: str = "*") -> Dict[str, str]:
    """Plain codon -> amino-acid dict in the formats the analysis modules use.

    ``rna=True`` writes U instead of T; ``stop_symbol`` replaces "*" (several
    modules use "STOP"). Returns a new dict on every call.
    """
    table = get_genetic_code(table_id)
    result: Dict[str, str] = {}
    for codon, amino_acid in table.codons.items():
        key = codon.replace("T", "U") if rna else codon
        result[key] = stop_symbol if amino_acid == "*" else amino_acid
    return result
