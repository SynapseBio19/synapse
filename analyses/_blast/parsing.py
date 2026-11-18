"""BLAST XML (outfmt 5) parsing into the data model."""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from dataclasses import asdict
from typing import Any, Dict, Iterable, List, Optional, Tuple

from analyses._blast.config import DATABASES, PROGRAMS
from analyses._blast.models import (
    BlastAnalysisError,
    BlastHit,
    BlastHSP,
    BlastParameters,
    BlastQuery,
    BlastRun,
    BlastRunMetadata,
)


def _local_tag(element: ET.Element) -> str:
    """Return XML tag without namespace."""
    return element.tag.rsplit("}", 1)[-1]


def _child_text(element: ET.Element, tag: str, default: str = "") -> str:
    for child in element:
        if _local_tag(child) == tag:
            return (child.text or default).strip()
    return default


def _descendants(element: ET.Element, tag: str) -> Iterable[ET.Element]:
    for node in element.iter():
        if _local_tag(node) == tag:
            yield node


def _to_int(value: str, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _to_float(value: str, default: Optional[float] = None) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


# ============================================================================
# BLAST XML PARSER
# ============================================================================

def _evalue_or_infinity(text: str) -> float:
    """The E-value in ``text``; infinity if it is missing or not a number."""
    value = _to_float(text)
    return math.inf if value is None or math.isnan(value) else value


def _parse_hsp(hsp_element: ET.Element) -> BlastHSP:
    return BlastHSP(
        number=_to_int(_child_text(hsp_element, "Hsp_num"), 1) or 1,
        bit_score=_to_float(_child_text(hsp_element, "Hsp_bit-score"), 0.0) or 0.0,
        raw_score=_to_float(_child_text(hsp_element, "Hsp_score")),
        # BLAST reports E = 0 for its strongest matches. A missing or unreadable
        # value is infinite; 0 must stay 0 (it became infinity with "or" until
        # 2026-11-18, which ranked perfect hits last as NOT_SIGNIFICANT).
        evalue=_evalue_or_infinity(_child_text(hsp_element, "Hsp_evalue")),
        query_start=_to_int(_child_text(hsp_element, "Hsp_query-from"), 0) or 0,
        query_end=_to_int(_child_text(hsp_element, "Hsp_query-to"), 0) or 0,
        subject_start=_to_int(_child_text(hsp_element, "Hsp_hit-from"), 0) or 0,
        subject_end=_to_int(_child_text(hsp_element, "Hsp_hit-to"), 0) or 0,
        query_frame=_to_int(_child_text(hsp_element, "Hsp_query-frame")),
        subject_frame=_to_int(_child_text(hsp_element, "Hsp_hit-frame")),
        identity=_to_int(_child_text(hsp_element, "Hsp_identity"), 0) or 0,
        positive=_to_int(_child_text(hsp_element, "Hsp_positive"), 0) or 0,
        gaps=_to_int(_child_text(hsp_element, "Hsp_gaps"), 0) or 0,
        alignment_length=_to_int(_child_text(hsp_element, "Hsp_align-len"), 0) or 0,
        query_sequence=_child_text(hsp_element, "Hsp_qseq"),
        subject_sequence=_child_text(hsp_element, "Hsp_hseq"),
        midline=_child_text(hsp_element, "Hsp_midline"),
    )


def _infer_scientific_name_from_definition(definition: str) -> Optional[str]:
    """Infer a binomial organism name when BLAST XML omits taxonomy fields.

    Some NCBI BLAST XML responses provide a useful ``Hit_def`` but no
    ``Hit_sciname``/``Hit_taxid`` element. In that case Synapse may extract a
    conservative binomial-looking name from the definition for display. This
    is explicitly an inference from the returned definition, not a taxonomy
    lookup and not an authoritative NCBI taxonomy assignment.
    """
    text = re.sub(r"\s+", " ", (definition or "").strip())
    if not text:
        return None

    bracketed = re.findall(
        r"\[\s*([A-Z][A-Za-z-]{2,})\s+([a-z][A-Za-z-]{2,})\s*\]",
        text,
    )
    if bracketed:
        genus, species = bracketed[0]
        return f"{genus} {species}"

    candidates = re.findall(
        r"(?<![A-Za-z-])([A-Z][A-Za-z-]{2,})\s+([a-z][A-Za-z-]{2,})(?![A-Za-z-])",
        text,
    )
    if not candidates:
        return None

    excluded = {
        "Predicted",
        "Putative",
        "Hypothetical",
        "Uncharacterized",
        "Unnamed",
        "Novel",
        "Example",
        "Reference",
        "Sequence",
        "Protein",
        "Transcript",
    }
    for genus, species in candidates:
        if genus not in excluded:
            return f"{genus} {species}"

    return None


# The terminal dialogs of this module are in ui/blast_menu.py (P2-2).


def _parse_hit(hit_element: ET.Element, query_length: Optional[int] = None) -> BlastHit:
    hsp_container = next(
        (node for node in hit_element if _local_tag(node) == "Hit_hsps"),
        None,
    )
    hsps: List[BlastHSP] = []

    if hsp_container is not None:
        for hsp_node in hsp_container:
            if _local_tag(hsp_node) == "Hsp":
                hsps.append(_parse_hsp(hsp_node))

    definition = _child_text(hit_element, "Hit_def")
    scientific_names = []
    for tag in ("Hit_sciname", "Hit_scientific_name", "Hit_taxname"):
        value = _child_text(hit_element, tag)
        if value and value not in scientific_names:
            scientific_names.append(value)

    taxonomy_ids = []
    for tag in ("Hit_taxid", "Hit_taxid_list"):
        value = _child_text(hit_element, tag)
        if value:
            taxonomy_ids.extend(
                item.strip() for item in value.split(";") if item.strip()
            )

    # Legacy/compact NCBI XML may omit explicit taxonomy fields. Preserve the
    # real BLAST definition and add only a conservative, clearly documented
    # fallback label so the taxonomy overview remains useful.
    if not scientific_names:
        inferred_name = _infer_scientific_name_from_definition(definition)
        if inferred_name:
            scientific_names.append(inferred_name)

    return BlastHit(
        hit_id=_child_text(hit_element, "Hit_id"),
        accession=_child_text(hit_element, "Hit_accession"),
        definition=definition,
        length=_to_int(_child_text(hit_element, "Hit_len")),
        hsps=hsps,
        scientific_names=scientific_names,
        taxonomy_ids=taxonomy_ids,
        query_length=query_length,
    )

def _parse_parameters(root: ET.Element) -> Dict[str, Any]:
    parameters: Dict[str, Any] = {}
    for container in _descendants(root, "Parameters"):
        for child in container:
            tag = _local_tag(child)
            key = tag.replace("Parameters_", "")
            value = (child.text or "").strip()
            parameters[key] = value
        break
    return parameters


def _parse_database_statistics(root: ET.Element) -> Dict[str, Any]:
    stats: Dict[str, Any] = {}
    for container in _descendants(root, "Statistics"):
        for child in container:
            tag = _local_tag(child).replace("Statistics_", "")
            value = (child.text or "").strip()
            parsed: Any = value
            if re.fullmatch(r"-?\d+", value):
                parsed = int(value)
            elif re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)(?:[Ee][+-]?\d+)?", value):
                parsed = float(value)
            stats[tag] = parsed
        break
    return stats


def describe_run_from_xml(xml_text: str) -> Tuple[BlastQuery, BlastParameters]:
    """Rebuild a query and parameters from a BLAST XML document (P1-7).

    Used when a run is resumed from its request ID and the caller no longer
    holds the original objects. Only what the XML states is recovered: the
    query sequence is not part of BLAST output, so the query carries its
    reported length and an empty sequence rather than an invented one.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise BlastAnalysisError(f"Invalid BLAST XML: {exc}") from exc

    definition = _child_text(root, "BlastOutput_query-def") or "Resumed_Query"
    raw_length = _child_text(root, "BlastOutput_query-len") or "0"
    try:
        length = int(raw_length)
    except (TypeError, ValueError):
        length = 0
    program = (_child_text(root, "BlastOutput_program") or "blastn").lower().strip()
    database = _child_text(root, "BlastOutput_db") or "unknown"

    query = BlastQuery(
        name=definition.split()[0] if definition.split() else "Resumed_Query",
        sequence="",
        sequence_type=PROGRAMS.get(program, {}).get("query_type", "DNA"),
        source="Resumed NCBI request",
        region_start=1,
        region_end=max(1, length),
        feature_name=definition,
        notes=["Recovered from BLAST XML; the query sequence is not part of it."],
        reported_length=length,
    )
    if not length:
        query.notes.append(
            "The XML states no query length, so coverage cannot be computed "
            "for this resumed run."
        )
    return query, BlastParameters(program=program, database=database)


def parse_blast_xml(
    xml_text: str,
    *,
    query: BlastQuery,
    rid: str,
    parameters: BlastParameters,
    submitted_at: str,
    completed_at: Optional[str] = None,
    database_type_override: Optional[str] = None,
) -> BlastRun:
    """Parse NCBI BLAST XML into a Synapse BlastRun."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise BlastAnalysisError(f"Invalid BLAST XML: {exc}") from exc

    program = (
        _child_text(root, "BlastOutput_program")
        or parameters.program
    )
    version = _child_text(root, "BlastOutput_version") or None
    database_reported = _child_text(root, "BlastOutput_db") or None

    hits: List[BlastHit] = []

    # Legacy XML has Iteration -> Iteration_hits -> Hit.
    # Search generically to remain tolerant of wrapper/container differences.
    for hit_element in _descendants(root, "Hit"):
        hits.append(_parse_hit(hit_element, query_length=query.length))

    statistics = _parse_database_statistics(root)

    metadata = BlastRunMetadata(
        rid=rid,
        submitted_at_utc=submitted_at,
        completed_at_utc=completed_at,
        program=program,
        database=parameters.database,
        database_type=(
            database_type_override
            or DATABASES.get(parameters.database, {}).get("type", "unknown")
        ),
        query_name=query.name,
        query_length=query.length,
        query_source=query.source,
        parameters=asdict(parameters),
        ncbi_status="READY",
        ncbi_version=version,
        ncbi_database_reported=database_reported,
        reported_parameters=_parse_parameters(root),
        database_details={
            "reported_name": database_reported,
            "sequence_count": statistics.get("db-num"),
            "total_length": statistics.get("db-len"),
            "note": (
                "BLAST XML does not state the database release date; the "
                "search timestamp identifies the database state for remote runs."
            ),
        },
    )

    notes = [
        "BLAST reports local sequence similarity and statistical significance.",
        "A high-scoring hit does not by itself prove biological identity or function.",
        "Interpretation depends on database composition, search parameters and query context.",
    ]

    if not hits:
        notes.append("No BLAST hits were returned for the submitted search.")

    return BlastRun(
        metadata=metadata,
        query=query,
        hits=hits,
        database_statistics=statistics,
        interpretation_notes=notes,
    )
