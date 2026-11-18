"""Query construction and validation: sequence normalization, type inference,
program compatibility, FASTA queries and the ORF adapter."""

from __future__ import annotations

from core.alphabets import (
    AMBIGUOUS,
    classify_sequence,
    declared_type_from_molecule,
    INVALID,
    RNA_IUPAC,
)
from pathlib import Path
from typing import Any, List, Mapping, Optional, Sequence, Tuple

from analyses._blast.paths import _resolve_synapse_user_path
from analyses._blast.config import AMBIGUOUS_DNA, DNA_BASES, PROGRAMS, PROTEIN_AMINO_ACIDS
from analyses._blast.models import BlastAnalysisError, BlastQuery


def programs_for_query_type(query_type: str) -> List[str]:
    """BLAST programs that accept a query of this type ("DNA" or "PROTEIN")."""
    wanted = str(query_type or "").strip().upper()
    return [name for name, info in PROGRAMS.items() if info["query_type"] == wanted]


def check_program_query_type(program: str, query: "BlastQuery") -> None:
    """Raise BlastAnalysisError unless *program* accepts *query*'s type."""
    name = str(program or "").lower().strip()
    if name not in PROGRAMS:
        raise BlastAnalysisError(f"Unsupported BLAST program: {name}")
    expected = PROGRAMS[name]["query_type"]
    if expected != str(query.sequence_type).upper():
        raise BlastAnalysisError(
            f"{name} expects a {expected} query, "
            f"but the query is {query.sequence_type}."
        )


def normalize_sequence(sequence: str, sequence_type: str = "DNA") -> str:
    """Normalize a sequence while preserving biological alphabet semantics."""
    if not isinstance(sequence, str):
        raise BlastAnalysisError("Sequence must be a string.")

    cleaned = "".join(sequence.split()).upper()
    if not cleaned:
        raise BlastAnalysisError("Sequence is empty.")

    sequence_type = sequence_type.upper()
    if sequence_type == "DNA":
        allowed = DNA_BASES | AMBIGUOUS_DNA
        invalid = sorted(set(cleaned) - allowed)
        if invalid:
            raise BlastAnalysisError(
                f"Invalid DNA characters: {', '.join(invalid)}"
            )
    elif sequence_type == "RNA":
        invalid = sorted(set(cleaned) - RNA_IUPAC)
        if invalid:
            raise BlastAnalysisError(
                f"Invalid RNA characters: {', '.join(invalid)}"
            )
    elif sequence_type == "PROTEIN":
        invalid = sorted(set(cleaned) - PROTEIN_AMINO_ACIDS)
        if invalid:
            raise BlastAnalysisError(
                f"Invalid protein characters: {', '.join(invalid)}"
            )
    else:
        raise BlastAnalysisError("sequence_type must be DNA, RNA or PROTEIN.")

    return cleaned


def infer_sequence_type(sequence: str, declared: Optional[str] = None) -> str:
    """Return "DNA", "RNA" or "PROTEIN" using ``core.alphabets.classify_sequence``.

    Raises BlastAnalysisError when the type is invalid or cannot be told
    from the content (AMBIGUOUS); the caller must then pass it explicitly.
    *declared* may be a type or "NUCLEOTIDE" (e.g. from file metadata).
    """
    result = classify_sequence(sequence, declared)
    if result.evidence == INVALID:
        raise BlastAnalysisError(f"Cannot infer the query type: {result.reasons[0]}.")
    if result.evidence == AMBIGUOUS:
        raise BlastAnalysisError(
            "The query type cannot be determined from the sequence "
            f"({result.reasons[0]}); specify DNA or PROTEIN."
        )
    assert result.sequence_type is not None
    return result.sequence_type


def build_query(
    sequence: str,
    *,
    name: str = "Synapse_Query",
    sequence_type: Optional[str] = None,
    source: str = "Selected Sequence",
    region_start: int = 1,
    region_end: Optional[int] = None,
    feature_type: str = "sequence",
    feature_name: Optional[str] = None,
) -> BlastQuery:
    inferred = sequence_type.upper() if sequence_type else infer_sequence_type(sequence)
    normalized = normalize_sequence(sequence, inferred)
    notes: List[str] = []
    if inferred == "RNA":
        # BLAST nucleotide databases are searched in the DNA alphabet.
        normalized = normalized.replace("U", "T")
        inferred = "DNA"
        notes.append("RNA query searched as DNA: U read as T.")
    end = region_end if region_end is not None else region_start + len(normalized) - 1
    if region_start < 1 or end < region_start:
        raise BlastAnalysisError("Invalid query coordinates.")

    return BlastQuery(
        name=name.strip() or "Synapse_Query",
        sequence=normalized,
        sequence_type=inferred,
        source=source,
        region_start=region_start,
        region_end=end,
        feature_type=feature_type,
        feature_name=feature_name,
        notes=notes,
    )


def load_fasta_queries(
    fasta_path: str | Path,
    *,
    source: str = "FASTA File",
) -> List[BlastQuery]:
    """Load one or more BLAST queries from a FASTA file.

    The parser is intentionally small and dependency-free because Synapse only
    needs FASTA headers and sequence records for BLAST query preparation.
    Every record is validated through :func:`build_query`, so invalid DNA or
    protein characters are rejected before a BLAST search is started. An
    NCBI ``[moltype=...]`` header modifier declares the record's type; a
    record whose type cannot be told from its letters needs one.

    A multi-record FASTA is returned as multiple ``BlastQuery`` objects. The
    interactive menu lets the user explicitly choose which record to search,
    avoiding an implicit first-record-only behavior.
    """
    path = _resolve_synapse_user_path(fasta_path)
    if not path.is_file():
        raise BlastAnalysisError(f"FASTA file not found: {path}")

    records: List[Tuple[str, str]] = []
    current_header: Optional[str] = None
    current_sequence: List[str] = []

    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except OSError as exc:
        raise BlastAnalysisError(f"Could not read FASTA file: {exc}") from exc
    except UnicodeError as exc:
        raise BlastAnalysisError(f"FASTA file is not valid UTF-8: {exc}") from exc

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(current_sequence)))
            current_header = line[1:].strip() or f"FASTA_Record_{len(records) + 1}"
            current_sequence = []
        else:
            if current_header is None:
                raise BlastAnalysisError(
                    "Invalid FASTA file: sequence data appears before the first header."
                )
            current_sequence.append(line)

    if current_header is not None:
        records.append((current_header, "".join(current_sequence)))

    if not records:
        raise BlastAnalysisError("FASTA file contains no sequence records.")

    queries: List[BlastQuery] = []
    for index, (header, sequence) in enumerate(records, start=1):
        if not sequence:
            raise BlastAnalysisError(
                f"FASTA record {index} ('{header}') contains no sequence."
            )
        from fileio.loader import fasta_molecule_type

        try:
            query_type = infer_sequence_type(
                sequence, declared_type_from_molecule(fasta_molecule_type(header))
            )
        except BlastAnalysisError as exc:
            raise BlastAnalysisError(
                f"FASTA record {index} ('{header}'): {exc} "
                "Add [moltype=protein] or [moltype=DNA] to its header."
            ) from exc
        queries.append(
            build_query(
                sequence,
                name=header,
                sequence_type=query_type,
                source=source,
            )
        )

    return queries


# ============================================================================
# ORF / SYNAPSE INTEGRATION
# ============================================================================

def extract_orf_context(orf: Any) -> Optional[Tuple[str, int, str]]:
    """Adapt common Synapse ORF dict/object representations."""
    if orf is None:
        return None

    def get_value(keys: Sequence[str], default: Any = None) -> Any:
        if isinstance(orf, Mapping):
            for key in keys:
                if key in orf:
                    return orf[key]
        for key in keys:
            if hasattr(orf, key):
                return getattr(orf, key)
        return default

    sequence = get_value(("sequence", "dna", "orf_sequence", "seq"))
    if not sequence:
        return None

    try:
        sequence = normalize_sequence(str(sequence), "DNA")
    except BlastAnalysisError:
        return None

    start = get_value(
        ("start", "start_position", "genomic_start", "position"),
        1,
    )
    try:
        start = int(start)
    except (TypeError, ValueError):
        start = 1

    name = str(
        get_value(("name", "id", "orf_id", "label"), "Selected ORF")
    )

    return sequence, max(1, start), name


def build_orf_query(
    orf: Any,
    *,
    sequence_name: str = "Selected Sequence",
    parent_sequence: Optional[str] = None,
) -> BlastQuery:
    """Build a BLAST query from a selected ORF.

    The query is the ORF in its coding orientation. region_start/region_end
    are 1-based forward-strand positions on the parent sequence (Synapse ORF
    dicts store 0-based starts; audit P0-2). The strand is recorded in
    ``feature_name``.
    """
    from core.orf_context import resolve_orf_region

    region = resolve_orf_region(orf, parent_sequence)
    if region is None:
        raise BlastAnalysisError("The selected ORF could not be extracted.")

    try:
        sequence = normalize_sequence(region.coding_sequence, "DNA")
    except BlastAnalysisError as exc:
        raise BlastAnalysisError(f"The selected ORF could not be extracted: {exc}") from exc

    start = region.start_1based or 1
    feature_name = region.name
    if region.start_1based is not None:
        feature_name = f"{region.name} ({region.strand} strand)"

    return build_query(
        sequence,
        name=region.name,
        sequence_type="DNA",
        source="Selected ORF",
        region_start=start,
        region_end=start + len(sequence) - 1,
        feature_type="ORF",
        feature_name=feature_name,
    )
