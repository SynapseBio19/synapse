"""Split plan for analyses/protein.py (Phase 9)."""

SOURCE = "analyses/protein.py"
PACKAGE = "analyses._protein"
PACKAGE_DIR = "analyses/_protein"
PACKAGE_DOC = """Implementation of :mod:`analyses.protein` (split out in Phase 9).

Import from ``analyses.protein``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``constants``: alphabets, residue masses, hydropathy and flexibility
  scales, pKa values, the instability (DIWV) table
- ``sequence``: validation, resolved and positional views, sequence quality
- ``translation``: codon resolution (IUPAC), RNA translation
- ``composition``: amino-acid composition (``analyze_protein``), molecular weight, extinction
  coefficient, aromaticity, aliphatic and instability index
- ``charge``: net charge and isoelectric point
- ``profiles``: hydropathy, hydrophobic and transmembrane regions,
  secondary-structure propensity, flexibility, signal peptide
- ``report``: the combined analyses and the preview
"""

MODULES = {
    "constants": ("Amino-acid alphabets, residue masses, scales, pKa sets and the DIWV table.", [
        "STANDARD_AMINO_ACIDS", "AA_ORDER", "NONSTANDARD_RESIDUES", "AMBIGUOUS_RESIDUE_OPTIONS",
        "UNRESOLVED_RESIDUES", "STOP_SYMBOL", "RNA_IUPAC", "AMINO_ACID_RESIDUE_MASSES", "WATER_MASS",
        "HYDROPHOBIC_RESIDUES", "POLAR_RESIDUES", "POSITIVE_RESIDUES", "NEGATIVE_RESIDUES",
        "AROMATIC_RESIDUES", "ALIPHATIC_RESIDUES", "KYTE_DOOLITTLE", "FLEXIBILITY_SCALE", "POSITIVE_PKA",
        "NEGATIVE_PKA", "NTERM_PKA", "CTERM_PKA", "INSTABILITY_ORDER", "INSTABILITY_ROWS",
        "build_instability_weights", "INSTABILITY_WEIGHTS",
    ]),
    "sequence": ("Protein sequence validation, resolved and positional views, sequence quality.", [
        "_strip_terminal_stop", "validate_protein_sequence", "get_resolved_protein", "get_positional_protein",
        "analyze_protein_quality",
    ]),
    "translation": ("Codon resolution with IUPAC symbols and RNA translation.", [
        "expand_iupac_codon", "rna_codon_table", "resolve_codon", "translate_protein", "analyze_translation",
    ]),
    "composition": ("Composition, molecular weight, extinction coefficient and composition indices.", [
        "analyze_protein", "calculate_amino_acid_percentages", "calculate_protein_molecular_weight",
        "calculate_aliphatic_index",
        "calculate_aromaticity", "calculate_extinction_coefficient", "calculate_instability_index",
    ]),
    "charge": ("Net charge at a pH and the isoelectric point.", [
        "_get_pka_tables", "calculate_charge_at_ph", "calculate_isoelectric_point",
    ]),
    "profiles": ("Sequence profiles and region predictions (heuristics, see METHODS).", [
        "calculate_gravy", "calculate_hydropathy_profile", "find_hydrophobic_regions",
        "find_transmembrane_domains", "calculate_secondary_structure_propensity",
        "calculate_secondary_structure_fraction", "calculate_flexibility_profile", "predict_signal_peptide",
    ]),
    "report": ("The combined protein analyses and the sequence preview.", [
        "analyze_protein_properties", "analyze_complete_protein", "create_protein_preview",
    ]),
}

FACADE_NAMES = set()
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 9: the implementation moved to the internal package\n"
    "# ``analyses._protein`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface and re-exports every name it\n"
    "# defined before.\n"
)
