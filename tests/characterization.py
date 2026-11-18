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

To rewrite the snapshots (only with a reason in docs/CHANGES.md):
``SYNAPSE_UPDATE_CHARACTERIZATION=1 python -m unittest tests.test_characterization``;
a section name instead of ``1`` (for example ``mutation``) rewrites only
the file of that section. Later phases add sections (``mutation.json``,
recorded from 727208a before the mutation module was split).
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


# ---------------------------------------------------------------------------
# Mutation (added in Phase 8, recorded from 727208a before its split)
# ---------------------------------------------------------------------------

def _variants(reference: str) -> dict:
    """Variants of ``reference`` with known edits (positions are 0-based)."""
    sub = reference[:40] + ("A" if reference[40] != "A" else "C") + reference[41:]
    ins = reference[:90] + "GGT" + reference[90:]
    dele = reference[:150] + reference[156:]
    multi = sub[:200] + ("T" if sub[200] != "T" else "G") + sub[201:230] + sub[232:]
    iupac = reference[:60] + "N" + reference[61:]
    return {"substitution": sub, "insertion": ins, "deletion": dele, "multiple": multi, "iupac": iupac}


def mutation_scenarios() -> dict:
    from analyses import mutation as M

    reference = _puc19()[146:446]
    variants = _variants(reference)
    out: dict = {}
    results = {}
    for name, variant in variants.items():
        result = M.compare_sequences(reference, variant, reference_name="ref", variant_name=name)
        results[name] = result
        out[f"compare_{name}"] = _plain(result)
        out[f"report_{name}"] = M.build_mutation_report(result, include_alignment=True)
        out[f"advanced_report_{name}"] = M.advanced_mutation_report(result)
        out[f"json_{name}"] = M.result_to_json(result)
        out[f"validate_{name}"] = _plain(M.validate_mutation_result(result))
        out[f"map_{name}"] = M.mutation_map(result.mutations, len(reference))
        out[f"signature_{name}"] = _plain(M.mutation_signature(result.mutations))
        out[f"hgvs_{name}"] = [[M.mutation_hgvs_c(m, reference), M.mutation_hgvs_p(m), M.format_mutation(m, reference)]
                               for m in result.mutations]
        if result.mutations:
            out[f"allele_{name}"] = M.hgvs_c_allele(result.mutations, reference)
            out[f"context_{name}"] = _plain(M.extract_mutation_context(result.mutations[0], reference, variant))
        out[f"restriction_{name}"] = _plain(M.restriction_site_impact(
            reference, variant, {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "PstI": "CTGCAG"}))
        out[f"primers_{name}"] = _plain(M.primer_mutation_overlap(result.mutations, [("fwd", 30, 50), ("rev", 140, 160)]))
        out[f"regions_{name}"] = _plain(M.annotate_regions(result.mutations, {"a": (0, 100), "b": (100, 300)}))
    multiple = M.compare_multiple_variants(reference, list(variants.values()))
    out["compare_multiple"] = _plain(multiple)
    out["shared"] = _plain(M.shared_mutations(multiple))
    out["unique"] = _plain(M.unique_mutations(results["multiple"], [results["substitution"]]))
    out["spectrum"] = _plain(M.build_mutation_spectrum(results["multiple"].mutations))
    out["hotspots"] = _plain(M.find_mutation_hotspots(results["multiple"].mutations, len(reference), window_size=30,
                                                      minimum_mutations=1))
    out["statistics"] = _plain(M.calculate_mutation_statistics(results["multiple"].mutations,
                                                                reference_length=len(reference)))
    out["quality"] = _plain(M.sequence_quality(reference + "AAAAAAAAAAAAAAAAAAAA"))
    out["indel_left"] = _plain(M.normalize_indel_left("ACGTTTTACG", "ACGTTTACG", 6, "T", ""))
    out["alignment"] = _plain(M.global_align(reference[:80], variants["multiple"][:85]))
    out["similarity"] = _plain(M.sequence_similarity(reference, variants["deletion"]))
    out["detect"] = _plain(M.detect_mutations(reference, variants["insertion"]))
    out["coding_frame"] = _plain(M.assess_coding_frame("ATGAAACCCGGGTTTTAA"))
    out["proteins"] = _plain(M.compare_proteins("MKVLLAAGG", "MKVLPAAG"))
    out["advanced"] = _plain(M.create_advanced_analysis(results["multiple"]))
    out["self_test"] = _plain(M.run_mutation_self_test())
    return out


SECTIONS = {
    SNAPSHOT: {"primer": primer_scenarios, "blast": blast_scenarios},
    PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "mutation.json": {"mutation": mutation_scenarios},
}


def all_scenarios() -> dict:
    return {"primer": primer_scenarios(), "blast": blast_scenarios()}


def snapshot_for(section: str) -> Path:
    return next(path for path, sections in SECTIONS.items() if section in sections)


def write_snapshot(only: str = "") -> None:
    """Rewrite the snapshot files; ``only`` limits it to one section's file."""
    for path, sections in SECTIONS.items():
        if only and only not in sections:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {name: function() for name, function in sections.items()}
        path.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n", encoding="utf-8")
