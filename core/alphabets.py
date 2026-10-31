"""Sequence alphabets and sequence-type classification.

One source of truth for which symbols make a DNA, RNA or protein sequence,
and for deciding which of these a given sequence is.

Alphabets
---------
* Nucleotides: the IUPAC/IUB codes (NC-IUB, Cornish-Bowden 1985): A, C, G,
  T (DNA) or U (RNA), and the ambiguity codes R Y S W K M B D H V N.
* Protein: the 20 standard one-letter amino-acid codes, the ambiguity codes
  B (Asx), Z (Glx), J (Xle) and X (any), the encoded residues U
  (selenocysteine) and O (pyrrolysine), and ``*`` for a stop.

Gap symbols (``-``, ``.``) and digits belong to no sequence alphabet.
``normalize_sequence`` only removes whitespace and uppercases: it never
converts U to T and never drops a symbol.

Classification
--------------
Content alone can settle the type only in one direction. A symbol that
exists only in the protein alphabet (E F I J L O P Q X Z ``*``) rules out
DNA and RNA. The reverse is never true: every IUPAC nucleotide code is also
an amino-acid letter, so a nucleotide-looking sequence is at most
*probably* nucleotide. ``classify_sequence`` therefore reports one of:

* ``CONFIRMED`` – the type was declared (by the caller or by file metadata)
  and the sequence fits it, or the sequence contains protein-only symbols.
* ``PROBABLE`` – only A, C, G, T/U and N: read as nucleotide, though it is
  also a valid peptide (Ala, Cys, Gly, Thr/Sec, Asn).
* ``AMBIGUOUS`` – nucleotide ambiguity codes other than N, which are also
  amino-acid letters, or both T and U. No type is chosen; the caller must
  declare one. This is a categorical rule, not a percentage cut-off.
* ``INVALID`` – empty, symbols outside every alphabet, or symbols that do
  not fit the declared type.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Tuple

DNA = "DNA"
RNA = "RNA"
PROTEIN = "PROTEIN"
SEQUENCE_TYPES: Tuple[str, ...] = (DNA, RNA, PROTEIN)

# Declared-only value: "a nucleotide sequence", DNA or RNA by its letters.
# File formats such as GenBank write RNA molecules in the DNA alphabet.
NUCLEOTIDE = "NUCLEOTIDE"

CONFIRMED = "CONFIRMED"
PROBABLE = "PROBABLE"
AMBIGUOUS = "AMBIGUOUS"
INVALID = "INVALID"
CONFIDENCE_LEVELS: Tuple[str, ...] = (CONFIRMED, PROBABLE, AMBIGUOUS, INVALID)

NUCLEOTIDE_AMBIGUITY: FrozenSet[str] = frozenset("RYSWKMBDHVN")
DNA_UNAMBIGUOUS: FrozenSet[str] = frozenset("ACGT")
RNA_UNAMBIGUOUS: FrozenSet[str] = frozenset("ACGU")
DNA_IUPAC: FrozenSet[str] = DNA_UNAMBIGUOUS | NUCLEOTIDE_AMBIGUITY
RNA_IUPAC: FrozenSet[str] = RNA_UNAMBIGUOUS | NUCLEOTIDE_AMBIGUITY

PROTEIN_STANDARD: FrozenSet[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")
PROTEIN_AMBIGUITY: FrozenSet[str] = frozenset("BZJX")
PROTEIN_NONSTANDARD: FrozenSet[str] = frozenset("UO")
PROTEIN_STOP = "*"
PROTEIN_IUPAC: FrozenSet[str] = (
    PROTEIN_STANDARD | PROTEIN_AMBIGUITY | PROTEIN_NONSTANDARD | frozenset(PROTEIN_STOP)
)

ALPHABETS: Dict[str, FrozenSet[str]] = {DNA: DNA_IUPAC, RNA: RNA_IUPAC, PROTEIN: PROTEIN_IUPAC}
ALL_SYMBOLS: FrozenSet[str] = DNA_IUPAC | RNA_IUPAC | PROTEIN_IUPAC

# Symbols whose presence rules out DNA and RNA.
PROTEIN_ONLY: FrozenSet[str] = PROTEIN_IUPAC - DNA_IUPAC - RNA_IUPAC

# A sequence made only of these is read as nucleotide (PROBABLE).
_PLAIN_NUCLEOTIDE: Dict[str, FrozenSet[str]] = {DNA: frozenset("ACGTN"), RNA: frozenset("ACGUN")}

_GAP_SYMBOLS = frozenset("-.")
_MAX_LISTED_SYMBOLS = 10

_DECLARED_ALIASES = {
    "DNA": DNA,
    "RNA": RNA,
    "PROTEIN": PROTEIN,
    "AA": PROTEIN,
    "NUCLEOTIDE": NUCLEOTIDE,
}


class SequenceTypeError(ValueError):
    """The sequence's type is invalid, undecided, or not accepted here."""

    def __init__(self, message: str, classification: "SequenceClassification") -> None:
        super().__init__(message)
        self.classification = classification


@dataclass(frozen=True)
class SequenceClassification:
    """Outcome of ``classify_sequence``.

    ``sequence_type`` is None for AMBIGUOUS and INVALID results;
    ``candidates`` lists every type whose alphabet the sequence fits.
    ``invalid_symbols`` holds (1-based position, symbol) pairs.
    """

    sequence_type: Optional[str]
    confidence: str
    candidates: Tuple[str, ...] = ()
    reasons: Tuple[str, ...] = ()
    declared: Optional[str] = None
    invalid_symbols: Tuple[Tuple[int, str], ...] = field(default=(), repr=False)

    @property
    def usable(self) -> bool:
        """True when a type was determined (CONFIRMED or PROBABLE)."""
        return self.confidence in (CONFIRMED, PROBABLE)

    def describe(self) -> str:
        label = self.sequence_type or "undetermined"
        text = f"{label} ({self.confidence.lower()})"
        return f"{text}: {'; '.join(self.reasons)}" if self.reasons else text

    def as_dict(self) -> Dict[str, object]:
        return {
            "sequence_type": self.sequence_type,
            "confidence": self.confidence,
            "candidates": list(self.candidates),
            "declared": self.declared,
            "reasons": list(self.reasons),
            "invalid_symbol_count": len(self.invalid_symbols),
        }


def normalize_sequence(sequence: str) -> str:
    """Remove whitespace and uppercase. Never converts or drops symbols."""
    return "".join(str(sequence or "").split()).upper()


def normalize_declared_type(value: Optional[str]) -> Optional[str]:
    """Map a declared type ("dna", "rna", "protein", "aa", "nucleotide") to its constant."""
    if value is None:
        return None
    key = str(value).strip().upper()
    if not key:
        return None
    if key not in _DECLARED_ALIASES:
        raise ValueError(
            f"Unknown sequence type {value!r}; expected DNA, RNA, PROTEIN or NUCLEOTIDE."
        )
    return _DECLARED_ALIASES[key]


def declared_type_from_molecule(label: Optional[str]) -> Optional[str]:
    """Declared type for a file's molecule label, or None if it names none.

    Protein labels ("aa", "protein", "peptide") give PROTEIN. Nucleotide
    labels (GenBank LOCUS "DNA", "mRNA", "ss-RNA", ..., NCBI FASTA
    ``moltype``/``molecule`` values such as "genomic", "dna", "rna") give
    NUCLEOTIDE: databases write RNA molecules in the DNA alphabet, so the
    letters decide between DNA and RNA.
    """
    text = str(label or "").strip().lower()
    if not text:
        return None
    if text in {"aa", "protein", "peptide", "prot"}:
        return PROTEIN
    if "dna" in text or "rna" in text or text in {"genomic", "nucleotide", "bp"}:
        return NUCLEOTIDE
    return None


def symbols_outside(sequence: str, alphabet: Iterable[str]) -> List[Tuple[int, str]]:
    """(1-based position, symbol) for every symbol not in *alphabet*."""
    allowed = frozenset(alphabet)
    return [(i, s) for i, s in enumerate(sequence, start=1) if s not in allowed]


def _listed(symbols: Iterable[str]) -> str:
    unique = sorted(set(symbols))
    shown = ", ".join(repr(s) for s in unique[:_MAX_LISTED_SYMBOLS])
    return shown + (", ..." if len(unique) > _MAX_LISTED_SYMBOLS else "")


def _invalid(reason: str, *, declared: Optional[str] = None,
             bad: Tuple[Tuple[int, str], ...] = (), candidates: Tuple[str, ...] = ()) -> SequenceClassification:
    return SequenceClassification(None, INVALID, candidates, (reason,), declared, bad)


def _classify_declared(sequence: str, declared: str, symbols: FrozenSet[str],
                       candidates: Tuple[str, ...]) -> SequenceClassification:
    if declared == NUCLEOTIDE:
        if symbols <= DNA_IUPAC:
            resolved = DNA
        elif symbols <= RNA_IUPAC:
            resolved = RNA
        else:
            bad = tuple(symbols_outside(sequence, DNA_IUPAC | RNA_IUPAC))
            if "T" in symbols and "U" in symbols:
                reason = "declared nucleotide, but contains both T and U"
            else:
                reason = f"declared nucleotide, but contains {_listed(s for _, s in bad)}"
            return _invalid(reason, declared=declared, bad=bad, candidates=candidates)
        letter = "U" if resolved == RNA else "T" if "T" in symbols else None
        detail = f"; written with {letter}" if letter else "; contains neither T nor U, read as DNA"
        return SequenceClassification(resolved, CONFIRMED, candidates,
                                      (f"declared nucleotide{detail}",), declared)

    bad = tuple(symbols_outside(sequence, ALPHABETS[declared]))
    if bad:
        reason = f"declared {declared}, but contains {_listed(s for _, s in bad)}"
        return _invalid(reason, declared=declared, bad=bad, candidates=candidates)
    return SequenceClassification(declared, CONFIRMED, candidates, (f"declared {declared}",), declared)


def classify_sequence(sequence: str, declared: Optional[str] = None) -> SequenceClassification:
    """Classify *sequence* as DNA, RNA or protein (see the module docstring).

    *declared* is a type stated by the caller or by file metadata (DNA, RNA,
    PROTEIN, or NUCLEOTIDE for "DNA or RNA"). The sequence is normalized
    with ``normalize_sequence`` first.
    """
    declared_type = normalize_declared_type(declared)
    seq = normalize_sequence(sequence)
    if not seq:
        return _invalid("empty sequence", declared=declared_type)

    symbols = frozenset(seq)
    foreign = tuple(symbols_outside(seq, ALL_SYMBOLS))
    if foreign:
        reason = f"contains {_listed(s for _, s in foreign)}: symbols of no sequence alphabet"
        if symbols & _GAP_SYMBOLS:
            reason += " (alignment gaps are not sequence symbols)"
        return _invalid(reason, declared=declared_type, bad=foreign)

    candidates = tuple(t for t in SEQUENCE_TYPES if symbols <= ALPHABETS[t])
    if declared_type is not None:
        return _classify_declared(seq, declared_type, symbols, candidates)

    protein_only = symbols & PROTEIN_ONLY
    if protein_only:
        return SequenceClassification(
            PROTEIN, CONFIRMED, candidates,
            (f"contains {_listed(protein_only)}: symbols that are not nucleotide codes",),
        )

    has_t, has_u = "T" in symbols, "U" in symbols
    if has_t and has_u:
        return SequenceClassification(
            None, AMBIGUOUS, candidates,
            ("contains both T and U, so it is neither DNA nor RNA; it is a valid peptide "
             "only if U is selenocysteine; declare the sequence type",),
        )

    nucleotide = RNA if has_u else DNA
    if symbols <= _PLAIN_NUCLEOTIDE[nucleotide]:
        bases = "A, C, G, U and N" if nucleotide == RNA else "A, C, G, T and N"
        reasons = [f"only {bases}; also a valid peptide, so the type is probable, not confirmed"]
        if not has_t and not has_u:
            reasons.append("contains neither T nor U, read as DNA")
        return SequenceClassification(nucleotide, PROBABLE, candidates, tuple(reasons))

    codes = symbols & (NUCLEOTIDE_AMBIGUITY - {"N"})
    return SequenceClassification(
        None, AMBIGUOUS, candidates,
        (f"contains {_listed(codes)}: nucleotide ambiguity codes that are also "
         "amino-acid letters; declare the sequence type",),
    )


def require_sequence_type(sequence: str, accepted: Iterable[str], *,
                          declared: Optional[str] = None,
                          context: str = "This analysis") -> SequenceClassification:
    """Classify *sequence* and raise SequenceTypeError unless its type is in *accepted*.

    INVALID and AMBIGUOUS results always raise; an AMBIGUOUS sequence is
    accepted once its type is declared.
    """
    accepted_types = tuple(accepted)
    result = classify_sequence(sequence, declared)
    wanted = " or ".join(accepted_types)
    if result.confidence == INVALID:
        raise SequenceTypeError(f"{context} needs a {wanted} sequence: {result.reasons[0]}.", result)
    if result.confidence == AMBIGUOUS:
        options = ", ".join(result.candidates) or "none"
        raise SequenceTypeError(
            f"{context} needs a {wanted} sequence, and this sequence's type cannot be "
            f"determined from its content (consistent with: {options}): "
            f"{result.reasons[0]}.", result)
    if result.sequence_type not in accepted_types:
        raise SequenceTypeError(
            f"{context} needs a {wanted} sequence, but this is {result.describe()}.", result)
    return result
