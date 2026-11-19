
# Phase 9: the implementation moved to the internal package
# ``analyses._orf`` (one module per responsibility, see its __init__).
# This module remains the public interface and re-exports every name it
# defined before.

from __future__ import annotations

from analyses._orf.codes import (
    ALL_START_CODONS,
    ALTERNATIVE_START_CODONS,
    analyze_codon_identity,
    classify_start_codon,
    codon_table_for,
    DEFAULT_MIN_ORF_LENGTH,
    DEFAULT_RANKING,
    DEFAULT_REPORT_NESTED,
    DEFAULT_START_CODON_POLICY,
    DNA_IUPAC,
    expand_iupac_codon,
    GENETIC_CODE,
    _genetic_code_name,
    initiation_codons,
    _is_alternative_start,
    is_initiation_codon,
    LEGACY_ORF_SETTINGS,
    ORF_RANKING_MODES,
    resolve_start_codons,
    STANDARD_START_CODONS,
    STANDARD_STOP_CODONS,
    START_CODON_POLICIES,
    start_codons_for,
    _stop_codons_of,
)
from analyses._orf.features import (
    analyze_codon_usage,
    analyze_kozak_context,
    analyze_orf_ambiguity,
    analyze_shine_dalgarno,
    analyze_start_context,
    calculate_cds_quality,
    calculate_codon_diversity,
    calculate_orf_confidence,
    calculate_orf_heuristic_score,
)
from analyses._orf.finder import (
    find_orfs,
    _five_prime_partial_orfs,
    _keep_longest_per_stop,
    MISSING_START_STATUSES,
    orf_end_on_sequence,
    orf_has_start_codon,
    _orf_sort_key,
    _origin_spanning_orfs,
)
from analyses._orf.summary import (
    analyze_orf_relationships,
    analyze_orfs,
    calculate_frame_statistics,
    calculate_orf_density,
    compare_orf_candidates,
    create_orf_preview,
    _orf_genomic_interval,
    rank_orf_candidates,
    summarize_orfs,
)
from analyses._orf.selection import (
    analyze_orfs_enhanced,
    analyze_selected_orf,
    calculate_orf_feature_summary,
    translate_orf_sequence,
)


# ============================================================
# ORF COMPARISON REPORT
# ============================================================


# ============================================================
# ORF UI / SELECTION
# ============================================================


# Backward-compatible UI aliases
