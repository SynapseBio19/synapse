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


# ---------------------------------------------------------------------------
# Cloning (added in Phase 8, recorded from 8e1f6d2 before its split)
# ---------------------------------------------------------------------------

def cloning_scenarios() -> dict:
    from analyses import cloning as C
    from tests._support import deterministic_dna

    vector = deterministic_dna(400) + "GAATTC" + deterministic_dna(300)[::-1] + "GGATCC" + deterministic_dna(200)
    insert = "GAATTC" + deterministic_dna(240) + "GGATCC"
    enzymes = ["EcoRI", "BamHI", "HindIII", "XhoI", "PstI", "SmaI", "KpnI", "XbaI"]
    out: dict = {}
    for topology in ("linear", "circular"):
        out[f"planner_{topology}"] = _plain(C.automatic_cloning_planner(
            vector, insert, enzyme_names=enzymes, top_n=5, vector_topology=topology))
        out[f"best_{topology}"] = _plain(C.find_best_cloning_strategies(
            vector, insert, enzyme_names=enzymes, top_n=5, vector_topology=topology))
        out[f"restriction_map_{topology}"] = _plain(C.build_restriction_map(vector, enzymes, topology))
    directional = C.evaluate_cloning_strategy(vector, insert, "EcoRI", "BamHI", "EcoRI", "BamHI")
    out["evaluate_directional"] = _plain(directional)
    out["evaluate_mismatched"] = _plain(C.evaluate_cloning_strategy(vector, insert, "EcoRI", "BamHI", "PstI", "KpnI"))
    out["pair"] = _plain(C.analyze_enzyme_pair(vector, insert, "EcoRI", "BamHI", "EcoRI", "BamHI"))
    single = C.analyze_cloning(vector, insert, "EcoRI", "EcoRI")
    out["single_enzyme"] = _plain(single)
    out["summary"] = _plain(C.get_cloning_summary(single))
    out["junctions"] = _plain(C.get_cloning_junctions(single))
    out["reasons"] = _plain(C.get_design_reasons(single))
    out["rank"] = _plain(C.rank_cloning_strategy(single))
    out["compare"] = _plain(C.compare_cloning_strategies(single, directional))
    out["suitability"] = {name: _plain(C.assess_enzyme_suitability(vector, insert, name)) for name in enzymes}
    out["compatible"] = _plain(C.find_compatible_enzymes(vector, insert, enzymes))
    for mode in ("auto", "advanced"):
        try:
            out[f"run_{mode}"] = _plain(C.run_cloning_analysis(vector, insert, mode=mode, enzyme_names=enzymes))
        except Exception as error:  # record the documented rejection as well
            out[f"run_{mode}"] = f"{type(error).__name__}: {error}"
    out["best_strategy"] = _plain(C.get_best_cloning_strategy(vector, insert, enzymes))
    out["overhangs"] = {name: _plain(C.get_overhang(name)) for name in enzymes}
    left = C.make_site(C.analyze_enzyme(vector, "EcoRI"), 400)
    right = C.make_site(C.analyze_enzyme(insert, "EcoRI"), 0)
    out["end_compatibility"] = _plain(C.analyze_end_compatibility(left, right))
    out["junction"] = _plain(C.analyze_junction(left, right, "left"))
    out["validation"] = [C.validate_dna_sequence("acgt nacg"), _plain(C.analyze_enzyme("GAATTCGAATTC", "EcoRI", "circular"))]
    out["self_test"] = _plain(C.run_cloning_self_test())
    return out


# ---------------------------------------------------------------------------
# Protein (added in Phase 9, recorded from ee8203d before its split)
# ---------------------------------------------------------------------------

PROTEINS = (
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAQSDILRRTQ",
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
    "MKXBZJACDEF*",
    "MLLAVLYCLLWSFQTSAGHFPRACVSSKNLMEKECCPPWSGDRSPCGQLSGRGSCQNILLSNAPLGPQFPFTGVDDRESWPSVFYNRTCQCSGNFMGFNCGNCKFGFWGPNCTERRLLVRRNIFDLSAPEKDKFFAYLTLAKHTISSDYVIPIGTYGQMKNGSTPMFNDINIYDLFVWMHYYVSMDALLGGSEIWRDIDFAHEAPAFLPWHRLFLLRWEQEIQKLTGDENFTIPYWDWRDAEKCDICTDEYMGGQHPTNPNLLSPASFFSSWQIVCSRLEEYNSHQSLCNGTPEGPLRRNPGNHDKSRTPRLPSSADVEFCLSLTQYESGSMDKAANFSFRNTLEGFASPLTGIADASQSSMHNALHIYMNGTMSQVQGSANDPIFLLHHAFVDSIFEQWLRRHRPLQEVYPEANAPIGHNRESYMVPFIPLYRNGDFFISSKDLGYDYSYLQDSDPDSFQDYIKSYLEQASRIWSWLLGAAMVGAVLTALLAGLVSLLCRHKRKQLPEEKQPLLMEKEDYHSLYQSHL",
)
RNAS = ("AUGGCUUCUAAAUAA", "AUGNNNGCUUAGGCU", "GUGAAACCCGGGUUUUGA", "AUGRYSWKMBDHVNUAA")


def protein_scenarios() -> dict:
    from analyses import protein as P

    out: dict = {}
    for index, sequence in enumerate(PROTEINS):
        out[f"complete_{index}"] = _plain(P.analyze_complete_protein(sequence))
        out[f"preview_{index}"] = P.create_protein_preview(sequence)
        out[f"quality_{index}"] = _plain(P.analyze_protein_quality(sequence))
        for name in ("calculate_amino_acid_percentages", "calculate_protein_molecular_weight",
                     "calculate_isoelectric_point", "calculate_gravy", "calculate_hydropathy_profile",
                     "find_hydrophobic_regions", "find_transmembrane_domains", "calculate_aliphatic_index",
                     "calculate_aromaticity", "calculate_extinction_coefficient", "calculate_instability_index",
                     "calculate_secondary_structure_propensity", "calculate_secondary_structure_fraction",
                     "calculate_flexibility_profile", "predict_signal_peptide", "analyze_protein_properties"):
            try:
                out[f"{name}_{index}"] = _plain(getattr(P, name)(sequence))
            except Exception as error:  # the documented rejection is part of the behavior
                out[f"{name}_{index}"] = f"{type(error).__name__}: {error}"
        out[f"charge_{index}"] = [repr(P.calculate_charge_at_ph(sequence, pH)) for pH in (3.0, 7.0, 7.4, 11.0)]
    for index, rna in enumerate(RNAS):
        for code in (1, 2, 11):
            out[f"translation_{index}_{code}"] = _plain(P.analyze_translation(rna, initiator_start=True,
                                                                              genetic_code=code))
        out[f"translate_{index}"] = _plain(P.translate_protein(rna))
        out[f"complete_rna_{index}"] = _plain(P.analyze_complete_protein(P.translate_protein(rna), rna=rna,
                                                                          initiator_start=True))
    out["codons"] = {codon: _plain(P.resolve_codon(codon)) for codon in ("AUG", "GCN", "RAY", "NNN", "UAR", "UGA")}
    out["iupac"] = _plain(P.expand_iupac_codon("RYN"))
    out["instability_weights_size"] = len(P.build_instability_weights())
    return out


# ---------------------------------------------------------------------------
# Variant analysis (added in Phase 9, recorded from 6e52ca7 before its split)
# ---------------------------------------------------------------------------

def _fasta_contigs(path: Path) -> dict:
    contigs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            contigs[name] = []
        elif name:
            contigs[name].append(line.strip())
    return {key: "".join(value) for key, value in contigs.items()}


def variant_scenarios() -> dict:
    from analyses import variant_analysis as V

    out: dict = {}
    norm = PROJECT_ROOT / "tests" / "fixtures" / "bcftools_norm"
    texts = {
        "test_variants": (PROJECT_ROOT / "analyses" / "test_variants.vcf").read_text(),
        "test_export": (PROJECT_ROOT / "analyses" / "test_export.vcf").read_text(),
        "bcftools_norm": (norm / "norm.vcf").read_text(),
        "bcftools_split": (norm / "norm.split.vcf").read_text(),
    }
    contigs = _fasta_contigs(norm / "norm.fa")
    for name, text in texts.items():
        result = V.analyze_vcf_text(text, name=name)
        out[f"analyze_{name}"] = _plain(result)
        out[f"report_{name}"] = V.build_variant_report(result)
        variants, metadata, samples = V.parse_vcf(text)
        out[f"vcf_{name}"] = V.variants_to_vcf(variants, metadata=metadata, sample_names=samples)
        out[f"tsv_{name}"] = V.variants_to_tsv(variants)
        out[f"statistics_{name}"] = _plain(V.calculate_variant_statistics(variants))
        out[f"cohort_{name}"] = _plain(V.calculate_cohort_summary(variants))
        out[f"priority_{name}"] = _plain(V.prioritize_variants(variants))
        for mode in ("exact", "chr_prefix"):
            try:
                out[f"normalize_{name}_{mode}"] = _plain(V.normalize_variants(
                    variants, contig_references=contigs, contig_matching=mode))
            except Exception as error:
                out[f"normalize_{name}_{mode}"] = f"{type(error).__name__}: {error}"
        filtered = V.filter_variants(variants, V.VariantFilter(minimum_quality=70.0))
        out[f"filter_{name}"] = _plain(filtered)
    left = V.parse_vcf(texts["test_variants"])[0]
    right = V.parse_vcf(texts["test_export"])[0]
    out["compare"] = _plain(V.compare_variant_sets(left, right))
    out["interval"] = _plain(V.variants_in_interval(left, "chr1", 150, 350))
    out["density"] = _plain(V.interval_variant_density(left, "chr1", 1, 1000))
    out["features"] = _plain(V.check_features_against_variants(
        [{"name": "geneA", "start": 90, "end": 250, "chromosome": "chr1"}], left))
    out["self_test"] = _plain(V.run_variant_analysis_self_test())
    return out


# ---------------------------------------------------------------------------
# Sequence overview (added in Phase 9, recorded from 94c2c0d before its split)
# ---------------------------------------------------------------------------

def sequence_scenarios() -> dict:
    from analyses import sequence as S

    puc = _puc19()
    sequences = {
        "puc_900": puc[:900],
        "iupac": puc[100:400] + "NNRYKMSW" + puc[400:500],
        "cpg_rich": "CGCGGCGCGTACGCGCGGGCGCCCGCGATCG" * 12 + puc[:200],
        "homopolymers": "AAAAAAAATTTTTGGGGGGCCCC" * 5 + "ACGT" * 30,
        "short": "ACGTNACG",
    }
    out: dict = {}
    for name, dna in sequences.items():
        out[f"analyze_{name}"] = _plain(S.analyze_sequence(dna, profile_window=50, profile_step=5))
        out[f"summary_{name}"] = _plain(S.create_sequence_summary(dna))
        for function in ("count_nucleotides", "count_ambiguous_bases", "get_complement", "get_reverse_complement",
                         "analyze_sequence_quality", "calculate_gc_content", "calculate_at_content",
                         "calculate_dna_content_ranges", "calculate_gc_skew", "calculate_at_skew",
                         "calculate_sequence_complexity", "analyze_cpg", "calculate_dna_molecular_weight",
                         "calculate_average_dna_molecular_weight", "get_sequence_forms",
                         "calculate_sequence_statistics"):
            try:
                out[f"{function}_{name}"] = _plain(getattr(S, function)(dna))
            except Exception as error:
                out[f"{function}_{name}"] = f"{type(error).__name__}: {error}"
        out[f"homopolymers_{name}"] = _plain(S.find_homopolymers(dna, min_length=4))
        out[f"cpg_islands_{name}"] = _plain(S.find_cpg_islands(dna, window_size=100))
        out[f"gc_profile_{name}"] = _plain(S.analyze_gc_profile(dna, window_size=40, step=7))
        out[f"complexity_profile_{name}"] = _plain(S.analyze_complexity_profile(dna, window_size=40, step=7))
    out["self_test"] = _plain(S.run_sequence_self_test())
    return out


# ---------------------------------------------------------------------------
# ORF finder (added in Phase 9, recorded from b46ed3d before its split)
# ---------------------------------------------------------------------------

def orf_scenarios() -> dict:
    from analyses import orf as O

    puc = _puc19()
    out: dict = {}
    settings = {
        "default": {},
        "legacy": dict(O.LEGACY_ORF_SETTINGS),
        "table_starts_code11": {"start_codon_policy": "table_starts", "genetic_code": 11, "min_length": 150},
        "nested_code2": {"nested": True, "genetic_code": 2, "min_length": 90},
        "partial5": {"include_5prime_partial": True, "min_length": 120},
    }
    for name, options in settings.items():
        for topology in ("linear", "circular"):
            result = O.analyze_orfs(puc, topology=topology, **options)
            out[f"analyze_{name}_{topology}"] = _plain(result)
    enhanced = O.analyze_orfs_enhanced(puc, min_length=150, top_n=5, topology="circular")
    out["enhanced"] = _plain(enhanced)
    orfs = O.analyze_orfs(puc, min_length=150, topology="circular")["orfs"]
    out["relationships"] = _plain(O.analyze_orf_relationships(orfs))
    out["ranked"] = _plain(O.rank_orf_candidates(orfs))
    out["compare"] = _plain(O.compare_orf_candidates(orfs, top_n=3))
    out["feature_summary"] = _plain(O.calculate_orf_feature_summary(orfs))
    out["summary"] = _plain(O.summarize_orfs(orfs, sequence_length=len(puc)))
    for index, item in enumerate(orfs[:4]):
        out[f"selected_{index}"] = _plain(O.analyze_selected_orf(item))
        out[f"translation_{index}"] = _plain(O.translate_orf_sequence(item))
    iupac = puc[:300] + "ATGNNNRYKAAATAG" + puc[300:600]
    out["iupac"] = _plain(O.analyze_orfs(iupac, min_length=30))
    out["codon_tables"] = {str(code): _plain(O.codon_table_for(code)) for code in (1, 4, 28)}
    out["starts"] = {f"{code}_{policy}": sorted(O.start_codons_for(code, policy))
                     for code in (1, 11) for policy in O.START_CODON_POLICIES}
    return out


# ---------------------------------------------------------------------------
# Restriction analysis (added in Phase 9, recorded from 709f964 before its split)
# ---------------------------------------------------------------------------

def restriction_scenarios() -> dict:
    from analyses import restriction as R

    puc = _puc19()
    enzymes = ["EcoRI", "BamHI", "HindIII", "PstI", "SmaI", "KpnI", "XbaI", "SalI", "NdeI", "BsaI", "SacI"]
    enzymes = [name for name in enzymes if name in R.RESTRICTION_ENZYMES]
    iupac = puc[:1200] + "GAANTC" + puc[1200:1500]
    out: dict = {"enzymes": enzymes}
    for topology in ("linear", "circular"):
        out[f"analysis_{topology}"] = _plain(R.run_restriction_analysis(puc, enzymes, topology))
        out[f"map_{topology}"] = _plain(R.create_restriction_map(puc, enzymes, topology))
        out[f"summary_{topology}"] = _plain(R.summarize_restriction_map(puc, enzymes, topology))
        digest = R.perform_digest(puc, ["EcoRI", "PstI", "NdeI"], topology=topology)
        out[f"digest_{topology}"] = _plain(digest)
        out[f"gel_{topology}"] = _plain(R.get_gel_band_prediction(digest))
        out[f"boundaries_{topology}"] = _plain(R.get_physical_cut_boundaries(digest))
        out[f"unique_{topology}"] = _plain(R.find_unique_cutters(puc, enzymes, topology))
        out[f"non_{topology}"] = _plain(R.find_non_cutters(puc, enzymes, topology))
        out[f"pairs_{topology}"] = _plain(R.recommend_enzyme_pairs(puc, enzymes, 4, topology))
        out[f"double_{topology}"] = _plain(R.run_double_digest(puc, "EcoRI", "HindIII", topology))
        for name in enzymes:
            out[f"sites_{name}_{topology}"] = _plain(R.summarize_enzyme_sites(puc, name, topology))
        out[f"iupac_{topology}"] = _plain(R.perform_digest(iupac, ["EcoRI"], topology=topology,
                                                           include_possible=True))
    out["validate_database"] = _plain(R.validate_enzyme_database())
    out["compatible"] = _plain(R.get_cloning_compatible_sites(puc, puc[:900], "EcoRI"))
    out["self_test"] = _plain(R.run_restriction_self_test())
    return out


SECTIONS = {
    SNAPSHOT: {"primer": primer_scenarios, "blast": blast_scenarios},
    PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "mutation.json": {"mutation": mutation_scenarios},
    PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "cloning.json": {"cloning": cloning_scenarios},
    PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "protein.json": {"protein": protein_scenarios},
    PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "variant.json": {"variant": variant_scenarios},
    PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "sequence.json": {"sequence": sequence_scenarios},
    PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "orf.json": {"orf": orf_scenarios},
    PROJECT_ROOT / "tests" / "fixtures" / "characterization" / "restriction.json": {
        "restriction": restriction_scenarios},
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
