
# Phase 9: the implementation moved to the internal package
# ``analyses._protein`` (one module per responsibility, see its __init__).
# This module remains the public interface and re-exports every name it
# defined before.

from __future__ import annotations

from analyses._protein.constants import (
    AA_ORDER,
    ALIPHATIC_RESIDUES,
    AMBIGUOUS_RESIDUE_OPTIONS,
    AMINO_ACID_RESIDUE_MASSES,
    AROMATIC_RESIDUES,
    build_instability_weights,
    CTERM_PKA,
    FLEXIBILITY_SCALE,
    HYDROPHOBIC_RESIDUES,
    INSTABILITY_ORDER,
    INSTABILITY_ROWS,
    INSTABILITY_WEIGHTS,
    KYTE_DOOLITTLE,
    NEGATIVE_PKA,
    NEGATIVE_RESIDUES,
    NONSTANDARD_RESIDUES,
    NTERM_PKA,
    POLAR_RESIDUES,
    POSITIVE_PKA,
    POSITIVE_RESIDUES,
    RNA_IUPAC,
    STANDARD_AMINO_ACIDS,
    STOP_SYMBOL,
    UNRESOLVED_RESIDUES,
    WATER_MASS,
)
from analyses._protein.sequence import (
    analyze_protein_quality,
    get_positional_protein,
    get_resolved_protein,
    _strip_terminal_stop,
    validate_protein_sequence,
)
from analyses._protein.translation import (
    analyze_translation,
    expand_iupac_codon,
    resolve_codon,
    rna_codon_table,
    translate_protein,
)
from analyses._protein.composition import (
    analyze_protein,
    calculate_aliphatic_index,
    calculate_amino_acid_percentages,
    calculate_aromaticity,
    calculate_extinction_coefficient,
    calculate_instability_index,
    calculate_protein_molecular_weight,
)
from analyses._protein.charge import (
    calculate_charge_at_ph,
    calculate_isoelectric_point,
    _get_pka_tables,
)
from analyses._protein.profiles import (
    calculate_flexibility_profile,
    calculate_gravy,
    calculate_hydropathy_profile,
    calculate_secondary_structure_fraction,
    calculate_secondary_structure_propensity,
    find_hydrophobic_regions,
    find_transmembrane_domains,
    predict_signal_peptide,
)
from analyses._protein.report import (
    analyze_complete_protein,
    analyze_protein_properties,
    create_protein_preview,
)


# ============================================================
# REPORT
# ============================================================

# The terminal dialogs of this module are in ui/protein_menu.py (P2-2).


# ============================================================
# PROTEIN ANALYSIS UI - SYNAPSE PROFESSIONAL V5
# ============================================================
#
# The analysis engine above remains unchanged.
#
# Menu:
#   [1] Protein Overview
#   [2] Translation & Sequence
#   [3] Amino Acid Composition
#   [4] Physicochemical & Structural Properties
#   [5] Complete Protein Report
#   [0] Back
#
# ============================================================


# ============================================================
# 2 - TRANSLATION & SEQUENCE
# ============================================================


# ============================================================
# 3 - AMINO ACID COMPOSITION
# ============================================================


# ============================================================
# 4 - PHYSICOCHEMICAL & STRUCTURAL PROPERTIES
# ============================================================


# ============================================================
# 5 - COMPLETE PROTEIN REPORT
# ============================================================


# ============================================================
# MENU
# ============================================================


# ============================================================
# MAIN PROTEIN UI
# ============================================================


# ============================================================
# PIPELINE ADAPTER
# ============================================================
