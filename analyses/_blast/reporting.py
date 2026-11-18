"""Serialization, CSV/JSON/text exports, the text report and run comparison."""

from __future__ import annotations

import csv
import json
from core.textformat import field_line, major_rule, minor_rule
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict

from analyses._blast.paths import _resolve_synapse_user_path
from analyses._blast.models import (
    BlastComparison,
    BlastHit,
    BlastHSP,
    BlastRun,
    _evalue_text,
    rank_hits,
)
from analyses._blast.evidence import build_evidence_profile


# ============================================================================
# SERIALIZATION / EXPORT
# ============================================================================

def hsp_to_dict(hsp: BlastHSP) -> Dict[str, Any]:
    data = asdict(hsp)
    data["identity_percent"] = hsp.identity_percent
    data["positive_percent"] = hsp.positive_percent
    data["gap_percent"] = hsp.gap_percent
    data["query_coverage"] = hsp.query_coverage
    data["alignment_nongap_percent"] = hsp.alignment_nongap_percent
    return data


def hit_to_dict(hit: BlastHit) -> Dict[str, Any]:
    data = asdict(hit)
    data["best_bit_score"] = hit.best_bit_score
    data["best_evalue"] = hit.best_evalue
    data["best_identity_percent"] = hit.best_identity_percent
    data["best_coverage_percent"] = hit.best_coverage_percent
    data["best_hsp_coverage_percent"] = hit.best_hsp_coverage_percent
    data["total_alignment_length"] = hit.total_alignment_length
    data["hsps"] = [hsp_to_dict(hsp) for hsp in hit.hsps]
    return data


def run_to_dict(run: BlastRun) -> Dict[str, Any]:
    return {
        "metadata": asdict(run.metadata),
        "query": asdict(run.query),
        "hits": [hit_to_dict(hit) for hit in rank_hits(run.hits)],
        "database_statistics": run.database_statistics,
        "interpretation_notes": run.interpretation_notes,
        "evidence_profile": [
            asdict(profile) for profile in build_evidence_profile(run)
        ],
    }


def run_to_json(run: BlastRun, *, indent: int = 2) -> str:
    return json.dumps(
        run_to_dict(run),
        indent=indent,
        ensure_ascii=False,
    )


def export_json(run: BlastRun, path: str | Path) -> str:
    path = _resolve_synapse_user_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(run_to_json(run), encoding="utf-8")
    return str(path)


def export_csv(run: BlastRun, path: str | Path) -> str:
    path = _resolve_synapse_user_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "rank",
        "accession",
        "hit_id",
        "definition",
        "identity_percent",
        "query_coverage_percent",
        "evalue",
        "bit_score",
        "alignment_length",
        "evidence_score",
        "evidence_band",
    ]

    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        profiles = {
            profile.accession: profile
            for profile in build_evidence_profile(run)
        }

        for rank, hit in enumerate(rank_hits(run.hits), start=1):
            profile = profiles.get(hit.accession)
            writer.writerow(
                {
                    "rank": rank,
                    "accession": hit.accession,
                    "hit_id": hit.hit_id,
                    "definition": hit.definition,
                    "identity_percent": (
                        profile.identity_percent if profile else hit.best_identity_percent
                    ),
                    "query_coverage_percent": (
                        profile.query_coverage_percent if profile else hit.best_coverage_percent
                    ),
                    "evalue": hit.best_evalue,
                    "bit_score": hit.best_bit_score,
                    "alignment_length": (
                        profile.alignment_length
                        if profile
                        else hit.best_hsp.alignment_length if hit.best_hsp else 0
                    ),
                    "evidence_score": profile.evidence_score if profile else "",
                    "evidence_band": profile.evidence_band if profile else "",
                }
            )

    return str(path)


def build_text_report(run: BlastRun) -> str:
    lines = [
        major_rule(),
        "SYNAPSE BLAST REPORT",
        major_rule(),
        field_line('Query', f"{run.query.name}"),
        field_line('Query type', f"{run.query.sequence_type}"),
        field_line('Query length', f"{run.query.length}"),
        field_line('Query source', f"{run.query.source}"),
        field_line('Feature type', f"{run.query.feature_type}"),
        field_line('Program', f"{run.metadata.program}"),
        field_line('Database', f"{run.metadata.database}"),
        field_line('RID', f"{run.metadata.rid}"),
        field_line('Hits', f"{run.hit_count}"),
        field_line('Significant hits', f"{run.significant_hit_count}"),
        "",
        "TOP HITS",
        minor_rule(),
    ]

    for index, hit in enumerate(rank_hits(run.hits)[:20], start=1):
        lines.extend(
            [
                f"[{index}] {hit.accession or hit.hit_id}",
                f"    Definition: {hit.definition}",
                f"    Identity:   {hit.best_identity_percent:.2f}%",
                f"    Coverage:   {hit.best_coverage_percent:.2f}%",
                f"    E-value:    {_evalue_text(hit.best_evalue)}",
                f"    Bit score:  {hit.best_bit_score:.2f}",
                "",
            ]
        )

    lines.extend(
        [
            "SCIENTIFIC INTERPRETATION NOTES",
            minor_rule(),
            *[f"- {note}" for note in run.interpretation_notes],
            "",
            "RUN METADATA",
            minor_rule(),
            json.dumps(asdict(run.metadata), indent=2, ensure_ascii=False),
        ]
    )

    return "\r\n".join(lines)


def export_text(run: BlastRun, path: str | Path) -> str:
    path = _resolve_synapse_user_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # The report already ends its lines with CRLF; newline="" writes them as
    # they are (text mode turned each into CR CR LF before 2026-11-12, D24).
    path.write_text(build_text_report(run), encoding="utf-8", newline="")
    return str(path)


# ============================================================================
# RUN COMPARISON
# ============================================================================

def compare_blast_runs(run_a: BlastRun, run_b: BlastRun) -> BlastComparison:
    map_a = {
        hit.accession or hit.hit_id: hit
        for hit in run_a.hits
    }
    map_b = {
        hit.accession or hit.hit_id: hit
        for hit in run_b.hits
    }

    keys_a = set(map_a)
    keys_b = set(map_b)
    common = sorted(keys_a & keys_b)

    identity_differences: Dict[str, float] = {}
    coverage_differences: Dict[str, float] = {}

    for accession in common:
        identity_differences[accession] = round(
            map_a[accession].best_identity_percent
            - map_b[accession].best_identity_percent,
            4,
        )
        coverage_differences[accession] = round(
            map_a[accession].best_coverage_percent
            - map_b[accession].best_coverage_percent,
            4,
        )

    return BlastComparison(
        query_a=run_a.query.name,
        query_b=run_b.query.name,
        common_accessions=common,
        unique_to_a=sorted(keys_a - keys_b),
        unique_to_b=sorted(keys_b - keys_a),
        identity_differences=identity_differences,
        coverage_differences=coverage_differences,
    )
