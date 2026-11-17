"""Characterization scenarios for the BLAST and primer modules (Phase 7.1).

Phase 7.1 moves the implementation of ``analyses/blast.py`` and
``analyses/primer.py`` into internal packages. These scenarios record what
the modules return *before* that move: the snapshot in
``tests/fixtures/characterization/blast_primer.json`` was written from the
code at commit 9077f1e. ``tests/test_characterization.py`` compares the
current code with it, so a refactoring that changes a result, a report text
or an exported file fails.

The values are recorded with ``repr`` for floats, so the comparison is
exact. Paths are recorded relative to the project root.

To rewrite the snapshot (only with a reason in docs/CHANGES.md):
``SYNAPSE_UPDATE_CHARACTERIZATION=1 python -m unittest tests.test_characterization``
"""

from __future__ import annotations

import dataclasses
import json
import os
from pathlib import Path
import tempfile
from typing import Any
from unittest import mock

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SNAPSHOT = PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "blast_primer.json"


def _plain(value: Any) -> Any:
    """Return a JSON-safe, exact representation of ``value``."""
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {"__type__": type(value).__name__,
                **{f.name: _plain(getattr(value, f.name)) for f in dataclasses.fields(value)}}
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted(_plain(v) for v in value)
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, Path):
        return _relative(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    return repr(value)


def _relative(path: Path | str) -> str:
    path = Path(path)
    try:
        return "<root>/" + path.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return "<outside>/" + path.name


def _scrub(value: Any, work: str) -> Any:
    """Replace a temporary directory in recorded strings by ``<work>``."""
    if isinstance(value, dict):
        return {k: _scrub(v, work) for k, v in value.items()}
    if isinstance(value, list):
        return [_scrub(v, work) for v in value]
    if isinstance(value, str):
        for prefix in {str(Path(work).resolve()), str(work)}:
            value = value.replace(prefix, "<work>")
        return value.replace("\\", "/")
    return value


def _puc19() -> str:
    lines = (PROJECT_ROOT / "tests" / "fixtures" / "pUC19_L09137.2.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


# ---------------------------------------------------------------------------
# Primer
# ---------------------------------------------------------------------------

PRIMERS = ("GTAAAACGACGGCCAGT", "CAGGAAACAGCTATGAC", "GCGCGCATATATGCGCGC", "ATTTTTTTTTAAAAAAAT")


def primer_scenarios() -> dict:
    from analyses import primer as P

    sequence = _puc19()[:700]
    out: dict = {}

    report = P.design_primers(sequence, 150, 550)
    out["design_primers"] = _plain(report)
    out["design_primers_summary"] = P.get_primer_summary(report)
    out["summarize_primer_pair"] = [P.summarize_primer_pair(pair, rank) for rank, pair in
                                    enumerate(report.ranked_pairs[:3], start=1)]

    for method in ("nearest_neighbor", "wallace", "salt_adjusted"):
        config = P.ProfessionalPrimerConfig(tm_method=method)
        pro = P.design_primers_professional(sequence, 150, 550, config=config)
        out[f"professional_{method}"] = _plain(pro)
        out[f"professional_{method}_summary"] = _plain(P.get_professional_primer_summary(pro))
        out[f"professional_{method}_text"] = [P.summarize_professional_primer_pair(pair, rank)
                                              for rank, pair in enumerate(pro.ranked_pairs[:3], start=1)]
        out[f"professional_{method}_binding_sites"] = [_plain(P.primer_pair_binding_sites(pair.pair))
                                                       for pair in pro.ranked_pairs[:2]]

    high_salt = P.ProfessionalPrimerConfig(thermodynamics=P.HIGH_SALT_SCREENING_CONFIG)
    out["professional_high_salt"] = _plain(P.design_primers_professional(sequence, 150, 550, config=high_salt))

    out["restriction_site_choices"] = _plain(P.restriction_site_choices(_puc19(), "EcoRI"))
    choices = P.restriction_site_choices(_puc19(), "EcoRI")
    if choices:
        out["restriction_site_design"] = _plain(
            P.design_primers_for_restriction_site(_puc19(), "EcoRI", int(choices[0]["site_start"])))

    for primer in PRIMERS:
        out[f"analyze_primer_{primer}"] = _plain(P.analyze_primer(primer))
        out[f"analyze_primer_professional_{primer}"] = _plain(P.analyze_primer_professional(primer))
        out[f"compare_tm_methods_{primer}"] = _plain(P.compare_tm_methods(primer))
        out[f"structures_{primer}"] = _plain([P.analyze_hairpin(primer), P.analyze_self_dimer(primer),
                                              P.analyze_heterodimer(primer, PRIMERS[0])])
        out[f"nn_{primer}"] = _plain(P.calculate_tm_nearest_neighbor(primer, return_details=True))
        out[f"tm_advanced_{primer}"] = _plain(P.calculate_tm_advanced(primer))
        out[f"terminal_{primer}"] = _plain(P.calculate_terminal_stability(primer))
    for name, conditions in (("default", P.ThermodynamicConfig()),
                             ("low_salt", P.LOW_SALT_SCREENING_CONFIG),
                             ("santalucia_2004", P.ThermodynamicConfig(nn_parameter_set="santalucia_2004"))):
        out[f"thermodynamic_summary_{name}"] = _plain(P.get_thermodynamic_summary(conditions))
        out[f"duplex_{name}"] = _plain(P.calculate_primer_thermodynamics(PRIMERS[0], conditions=conditions))
        out[f"tm_advanced_{name}"] = _plain(P.calculate_tm_advanced(PRIMERS[1], thermodynamic_config=conditions))
    out["self_tests"] = {
        "run_primer_self_test": P.run_primer_self_test(),
        "run_thermodynamic_self_test": P.run_thermodynamic_self_test(),
        "run_primer_v3_self_test_final": P.run_primer_v3_self_test_final(),
        "run_primer_v4_screening_self_test": P.run_primer_v4_screening_self_test(),
    }
    return out


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def blast_scenarios() -> dict:
    from analyses import blast as B

    out: dict = {}
    query = B.build_query("ATGCGTACGTACGTACGTAC" + "A" * 80, name="characterization")
    run = B.parse_blast_xml(B.SAMPLE_XML, query=query, rid="CHAR", parameters=B.BlastParameters(),
                            submitted_at="2026-01-01T00:00:00+00:00")
    out["run"] = _plain(run)
    out["run_to_json"] = B.run_to_json(run)
    out["text_report"] = B.build_text_report(run)
    out["evidence_profile"] = _plain(B.build_evidence_profile(run))
    out["taxonomy"] = B.group_hits_by_taxonomy(run.hits)
    out["alignment"] = B.alignment_text(run.hits[0])
    out["comparison"] = _plain(B.compare_blast_runs(run, run))
    out["describe_run_from_xml"] = _plain(B.describe_run_from_xml(B.SAMPLE_XML))
    with tempfile.TemporaryDirectory() as work:
        for name, export in (("csv", B.export_csv), ("json", B.export_json), ("txt", B.export_text)):
            path = export(run, Path(work) / f"run.{name}")
            out[f"export_{name}"] = Path(path).read_text(encoding="utf-8")
        root = Path(work) / "db"
        manager = B.LocalBlastDatabaseManager(root)
        out["empty_manager_databases"] = _scrub(_plain(manager.list_databases()), work)
        out["empty_manager_scan"] = _scrub(_plain(manager.scan_database("core_nt")), work)
    out["normalize"] = [B.normalize_sequence("acgt nacg"), B.normalize_sequence("MKV*", "PROTEIN")]
    out["infer"] = [B.infer_sequence_type(s) for s in ("ACGTACGT", "MKVLLA", "ACGUACGU")]
    out["programs"] = {t: B.programs_for_query_type(t) for t in B.QUERY_TYPES}
    out["catalog"] = _plain(B.LOCAL_DATABASE_CATALOG)
    out["databases"] = _plain(B.DATABASES)
    out["default_root"] = _relative(B.DEFAULT_LOCAL_DB_ROOT)
    out["user_path_relative"] = _relative(B._resolve_synapse_user_path("relative/x.json")).replace(
        _relative(Path.cwd()) + "/", "<cwd>/")
    with mock.patch.dict(os.environ, {"PATH": ""}):
        manager_bins = B.BlastExecutableManager()._candidate_bins()
    out["executable_candidate_bins_in_project"] = sorted(
        _relative(p) for p in manager_bins if _relative(p).startswith("<root>/"))
    out["self_test"] = B.run_self_test(verbose=False)
    return out


def all_scenarios() -> dict:
    return {"primer": primer_scenarios(), "blast": blast_scenarios()}


def write_snapshot() -> None:
    SNAPSHOT.parent.mkdir(parents=True, exist_ok=True)
    SNAPSHOT.write_text(json.dumps(all_scenarios(), indent=1, sort_keys=True) + "\n", encoding="utf-8")
