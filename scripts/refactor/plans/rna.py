"""Split plan for analyses/rna.py (Phase 9)."""

SOURCE = "analyses/rna.py"
PACKAGE = "analyses._rna"
PACKAGE_DIR = "analyses/_rna"
PACKAGE_DOC = """Implementation of :mod:`analyses.rna` (split out in Phase 9).

Import from ``analyses.rna``; it re-exports every name defined here and is
the only public interface. The modules, in dependency order:

- ``bases``: RNA base tables, codon table, transcription, counts, complements
- ``composition``: percentages, GC/AU content, ranges, skews, profiles,
  ambiguity, quality, molecular weight, complexity, homopolymers
- ``coding``: codon usage, Kozak and Shine-Dalgarno context, start/stop
  landscape, frame translation
- ``structure``: the hairpin screen (a heuristic; the energy model is
  ``analyses.rna_structure``)
- ``summary``: the feature score, the summary, ``analyze_rna``, the preview
- ``pipeline``: the RNA source of the session pipeline
"""

MODULES = {
    "bases": ("RNA base tables, codon table, transcription, base counts and complements.", [
        "RNA_STANDARD_BASES", "RNA_IUPAC_BASES", "RNA_AMBIGUOUS_CODES", "RNA_COMPLEMENT",
        "RNA_NUCLEOTIDE_RESIDUE_MASSES", "RNA_CODON_TABLE", "rna_codon_table", "_orf_value", "orf_genetic_code",
        "transcribe_rna", "count_rna_bases", "get_rna_complement", "get_rna_reverse_complement",
    ]),
    "composition": ("Composition, content profiles, quality, molecular weight, complexity and homopolymers.", [
        "calculate_rna_base_percentages", "calculate_rna_gc_content", "calculate_rna_au_content",
        "calculate_rna_content_ranges", "calculate_rna_gc_skew", "calculate_rna_au_skew",
        "calculate_rna_gc_profile", "analyze_rna_ambiguity", "analyze_rna_quality",
        "calculate_rna_molecular_weight", "calculate_rna_complexity", "calculate_rna_complexity_profile",
        "find_rna_homopolymers", "find_rna_homopolymer_regions",
    ]),
    "coding": ("Codon usage, translation-initiation context, start/stop landscape and frame translation.", [
        "analyze_rna_codon_usage", "calculate_rna_codon_diversity", "analyze_rna_kozak_context",
        "analyze_rna_shine_dalgarno", "find_rna_start_stop_landscape", "translate_rna_frame",
        "analyze_rna_three_frames",
    ]),
    "structure": ("The RNA hairpin screen (heuristic) and its summary.", [
        "find_rna_hairpins", "summarize_rna_structure",
    ]),
    "summary": ("The feature score, the analysis summary, analyze_rna and the preview.", [
        "calculate_rna_feature_score", "summarize_rna_analysis", "analyze_rna", "create_rna_preview",
    ]),
    "pipeline": ("The RNA source of the session pipeline (selected ORF or sequence).", [
        "_clean_sequence", "_get_pipeline_context", "_get_pipeline_selected_orf",
        "resolve_rna_source", "_store_rna_in_pipeline", "analyze_rna_from_pipeline", "run_rna_pipeline",
    ]),
}

FACADE_NAMES = set()
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 9: the implementation moved to the internal package\n"
    "# ``analyses._rna`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface and re-exports every name it\n"
    "# defined before.\n"
)
