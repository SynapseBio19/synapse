"""Golden-master scenarios: deterministic snapshots of current module behavior.

Each scenario returns plain data. ``tests/test_golden_master.py`` compares the
data with ``tests/golden/<name>.json``. A difference means behavior changed.
Intentional, scientifically justified changes are accepted by regenerating
the snapshot (``SYNAPSE_UPDATE_GOLDEN=1``) and documenting the reason in
``docs/CHANGES.md``. Snapshots record what the software does, not what is
scientifically correct; known-answer tests live in other test modules.
"""

from __future__ import annotations

from dataclasses import asdict

from tests._support import PROJECT_ROOT, deterministic_dna, read_fasta_records

EXAMPLES = dict(read_fasta_records(PROJECT_ROOT / "data" / "examples" / "dna.fasta"))
REFERENCE = EXAMPLES["TEST_8_REFERENCE"]
RANDOM_600 = deterministic_dna(600, seed=2026)


def scenario_sequence():
    from analyses import sequence as S
    result = S.analyze_sequence(REFERENCE + "NNRY" + RANDOM_600[:300])
    keep = ("length", "nucleotide_counts", "gc_content", "at_content", "gc_content_range",
            "gc_skew", "at_skew", "sequence_complexity", "homopolymers", "cpg", "cpg_islands",
            "molecular_weight", "average_molecular_weight", "ambiguous_bases", "reverse_complement")
    return {k: result[k] for k in keep}


def _orf_view(orf):
    keys = ("rank", "strand", "frame", "start", "stop", "length", "protein_length",
            "start_codon", "stop_codon", "status", "sequence")
    view = {k: orf.get(k) for k in keys}
    view["heuristic_score"] = orf["heuristic_score"]["score"]
    view["cds_quality"] = orf["cds_quality"]["score"]
    return view


def scenario_orf():
    from analyses import orf as O
    out = {}
    for name in ("TEST_8_REFERENCE", "TEST_8_REVERSE", "TEST_8_SNP_NONSENSE"):
        out[name] = [_orf_view(o) for o in O.find_orfs(EXAMPLES[name])]
    random_orfs = O.analyze_orfs_enhanced(RANDOM_600)
    out["random600_count"] = random_orfs["candidate_count"]
    out["random600_top5"] = [_orf_view(o) for o in random_orfs["orfs"][:5]]
    out["translate_GTG"] = O.translate_orf_sequence("GTGAAATTTTAA")
    return out


def scenario_rna():
    from analyses import rna as R
    a = R.analyze_rna(REFERENCE)
    return {
        "rna": a["rna"], "gc_content": a["gc_content"], "au_content": a["au_content"],
        "hairpin_count": a["structure"]["hairpin_count"],
        "frames": [(f["frame"], f["protein"], f["complete"]) for f in a["translation_frames"]["frames"]],
        "feature_score": a["feature_score"],
    }


def scenario_protein():
    from analyses import protein as P
    ubq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    r = P.analyze_complete_protein(ubq)
    t = P.analyze_complete_protein("", rna=REFERENCE.replace("T", "U"))
    return {
        "ubiquitin": {k: r[k] for k in ("length", "molecular_weight", "pI", "gravy",
                                        "instability_index", "aliphatic_index",
                                        "aromaticity", "extinction_coefficient", "charges")},
        "translated_reference": t["protein"],
        "gug_start": P.analyze_translation("GUGAAAUUUUAA")["protein"],
    }


def scenario_restriction():
    from analyses import restriction as R
    plasmid_like = RANDOM_600[:200] + "GAATTC" + RANDOM_600[200:400] + "GGATCC" + RANDOM_600[400:] + "AAGCTT"
    out = {}
    for topo in ("linear", "circular"):
        d = R.perform_digest(plasmid_like, ["EcoRI", "BamHI", "HindIII", "ClaI", "TaqI"], topo)
        out[topo] = {"cut_positions": d["cut_positions"],
                     "fragments": [(f["start"], f["end"], f["length"]) for f in d["fragments"]]}
    out["clai_site"] = R.analyze_restriction_site("AAAAATCGATAAAAA", "ClaI", 4)
    out["n_run_digest_cuts"] = R.perform_digest("A" * 50 + "N" * 10 + "A" * 50,
                                                ["EcoRI", "BamHI", "HindIII"], "linear")["cut_count"]
    return out


def scenario_cloning():
    from analyses import cloning as C
    names = sorted(C.RESTRICTION_ENZYMES)
    matrix = {}
    for a in names:
        for b in names:
            sa = C.make_site(C.analyze_enzyme("AAAA" + C.get_recognition_sequence(a) + "AAAA", a),
                             C.analyze_enzyme("AAAA" + C.get_recognition_sequence(a) + "AAAA", a)["sites"][0])
            sb = C.make_site(C.analyze_enzyme("CCCC" + C.get_recognition_sequence(b) + "CCCC", b),
                             C.analyze_enzyme("CCCC" + C.get_recognition_sequence(b) + "CCCC", b)["sites"][0])
            matrix[f"{a}+{b}"] = C.ends_compatible(sa, sb)
    return {"compatibility": matrix,
            "overhangs": {n: C.get_overhang(n) for n in names}}


def scenario_mutation():
    from analyses import mutation as M
    out = {}
    for name, seq in EXAMPLES.items():
        if name == "TEST_8_REFERENCE":
            continue
        r = M.detect_mutations(REFERENCE, seq)
        c = M.analyze_coding_consequences(REFERENCE, seq)
        out[name] = {
            "mutations": [(m.position, m.reference, m.alternative, m.mutation_type) for m in r.mutations],
            "identity": r.identity_percent,
            "coding": {k: getattr(c, k) for k in ("synonymous_count", "missense_count", "nonsense_count",
                                                  "start_loss_count", "stop_loss_count", "frameshift",
                                                  "premature_stop_position")},
        }
    return out


def scenario_variant():
    from analyses import variant_analysis as V
    text = (PROJECT_ROOT / "analyses" / "test_variants.vcf").read_text(encoding="utf-8")
    variants, _, _ = V.parse_vcf(text)
    stats = V.calculate_variant_statistics(variants)
    normalized = V.left_normalize_variant(V.Variant("chr1", 7, "AG", "A"), reference_sequence="GGCAAAAGTC")
    return {"variants": [(v.chromosome, v.position, v.reference, v.alternative, v.variant_type, v.filter_status)
                         for v in variants],
            "statistics": asdict(stats),
            "normalize_7_AG_A": (normalized.position, normalized.reference, normalized.alternative)}


def scenario_crispr():
    from analyses import crispr as C
    result = C.run_crispr_analysis(RANDOM_600[:300])
    return {"count": result.candidate_count,
            "top10": [(c.guide, c.pam, c.strand, c.guide_start, c.guide_end, c.rank_score)
                      for c in result.candidates[:10]]}


def scenario_blast():
    from analyses import blast as B
    query = B.build_query("ATGCGTACGTACGTACGTAC" + "A" * 80, name="golden")
    run = B.parse_blast_xml(B.SAMPLE_XML, query=query, rid="GOLDEN", parameters=B.BlastParameters(),
                            submitted_at="2026-01-01T00:00:00+00:00")
    hit = run.hits[0]
    ev = B.calculate_evidence_score(hit)
    return {"coverage": hit.best_coverage_percent, "identity": hit.best_identity_percent,
            "evidence_score": ev.evidence_score, "band": ev.evidence_band}


def scenario_literature():
    from analyses import literature as L
    arts = [
        L.LiteratureArticle(title="CRISPR screening in yeast", pmid="1", doi="10.1/a", authors=["A. B. Smith"]),
        L.LiteratureArticle(title="CRISPR screening in yeast", doi="10.1/A", authors=["Anna B. Smith"]),
        L.LiteratureArticle(title="CRISPR screening in yeast", pmid="2"),
        L.LiteratureArticle(title="Untitled other study", abstract="crispr yeast"),
    ]
    merged = L.deduplicate_articles(arts)
    ranked = L.rank_articles(merged, "CRISPR yeast screening")
    return {"merged": [(a.title, a.pmid, a.doi, a.authors) for a in merged],
            "ranked_scores": [(a.title, a.relevance_score) for a in ranked]}


def scenario_primer():
    from analyses import primer as P
    template = deterministic_dna(400, seed=77)
    report = P.design_primers_professional(template, 0, len(template), config=P.ProfessionalPrimerConfig())
    pairs = []
    for item in report.ranked_pairs[:5]:
        p = item.pair
        pairs.append((p.forward.sequence, p.reverse.sequence, p.amplicon_length, p.category,
                      item.forward_thermodynamics.tm_C, item.reverse_thermodynamics.tm_C))
    tms = {s: P.calculate_tm(s) for s in ("AGCGGATAACAATTTCACACAGGA", "GTAAAACGACGGCCAGT",
                                          "TAATACGACTCACTATAGGG")}
    return {"pairs": pairs, "tm_default": tms}


SCENARIOS = {
    "sequence": scenario_sequence,
    "orf": scenario_orf,
    "rna": scenario_rna,
    "protein": scenario_protein,
    "restriction": scenario_restriction,
    "cloning": scenario_cloning,
    "mutation": scenario_mutation,
    "variant": scenario_variant,
    "crispr": scenario_crispr,
    "blast": scenario_blast,
    "literature": scenario_literature,
    "primer": scenario_primer,
}
