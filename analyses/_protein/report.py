"""The combined protein analyses and the sequence preview."""

from __future__ import annotations

from analyses._protein.constants import (
    AROMATIC_RESIDUES,
    HYDROPHOBIC_RESIDUES,
    NEGATIVE_RESIDUES,
    POLAR_RESIDUES,
    POSITIVE_RESIDUES,
)
from analyses._protein.sequence import analyze_protein_quality, get_resolved_protein
from analyses._protein.translation import analyze_translation
from analyses._protein.composition import (
    analyze_protein,
    calculate_aliphatic_index,
    calculate_amino_acid_percentages,
    calculate_aromaticity,
    calculate_extinction_coefficient,
    calculate_instability_index,
    calculate_protein_molecular_weight,
)
from analyses._protein.charge import calculate_charge_at_ph, calculate_isoelectric_point
from analyses._protein.profiles import (
    calculate_flexibility_profile,
    calculate_gravy,
    calculate_hydropathy_profile,
    calculate_secondary_structure_propensity,
    find_hydrophobic_regions,
    find_transmembrane_domains,
    predict_signal_peptide,
)


def analyze_protein_properties(
    protein
):

    resolved = get_resolved_protein(
        protein
    )

    length = len(resolved)

    if length == 0:

        return {
            "resolved_residues":
                0,

            "hydrophobic_percentage":
                0.0,

            "polar_percentage":
                0.0,

            "positive_percentage":
                0.0,

            "negative_percentage":
                0.0,

            "aromatic_percentage":
                0.0
        }

    return {
        "resolved_residues":
            length,

        "hydrophobic_percentage":
            sum(
                residue
                in HYDROPHOBIC_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100,

        "polar_percentage":
            sum(
                residue
                in POLAR_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100,

        "positive_percentage":
            sum(
                residue
                in POSITIVE_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100,

        "negative_percentage":
            sum(
                residue
                in NEGATIVE_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100,

        "aromatic_percentage":
            sum(
                residue
                in AROMATIC_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100
    }


# ============================================================
# COMPLETE PROTEIN ANALYSIS
# ============================================================


def analyze_complete_protein(
    protein,
    rna=None,
    initiator_start=False,
    genetic_code=1,
):

    if not isinstance(protein, str):
        protein = ""

    protein = "".join(
        protein.upper().split()
    )

    translation = None

    if rna is not None:

        translation = analyze_translation(
            rna,
            initiator_start=initiator_start,
            genetic_code=genetic_code,
        )

        protein = translation[
            "protein"
        ]

    quality = (
        analyze_protein_quality(
            protein
        )
    )

    composition = (
        analyze_protein(
            protein
        )
    )

    percentages = (
        calculate_amino_acid_percentages(
            protein
        )
    )

    properties = (
        analyze_protein_properties(
            protein
        )
    )

    molecular_weight = (
        calculate_protein_molecular_weight(
            protein,
            translation
        )
    )

    pI = (
        calculate_isoelectric_point(
            protein
        )
    )

    charges = {}

    for ph in (
        5.0,
        7.0,
        7.4,
        9.0
    ):

        charges[str(ph)] = (
            calculate_charge_at_ph(
                protein,
                ph
            )
        )

    gravy = calculate_gravy(
        protein
    )

    hydropathy = (
        calculate_hydropathy_profile(
            protein,
            window=9
        )
    )

    hydrophobic_regions = (
        find_hydrophobic_regions(
            protein,
            window=19,
            threshold=1.6
        )
    )

    tm_domains = (
        find_transmembrane_domains(
            protein,
            window=19,
            hydropathy_threshold=1.6,
            max_charged=2
        )
    )

    instability = (
        calculate_instability_index(
            protein
        )
    )

    aliphatic = (
        calculate_aliphatic_index(
            protein
        )
    )

    aromaticity = (
        calculate_aromaticity(
            protein
        )
    )

    extinction = (
        calculate_extinction_coefficient(
            protein
        )
    )

    secondary_structure = (
        calculate_secondary_structure_propensity(
            protein
        )
    )

    flexibility = (
        calculate_flexibility_profile(
            protein,
            window=9
        )
    )

    signal_peptide = (
        predict_signal_peptide(
            protein
        )
    )

    return {
        "protein":
            protein,

        "length":
            len(protein),

        "quality":
            quality,

        "composition":
            composition,

        "percentages":
            percentages,

        "properties":
            properties,

        "molecular_weight":
            molecular_weight,

        "pI":
            pI,

        "charges":
            charges,

        "gravy":
            gravy,

        "hydropathy":
            hydropathy,

        "hydrophobic_regions":
            hydrophobic_regions,

        "transmembrane_domains":
            tm_domains,

        "instability_index":
            instability,

        "aliphatic_index":
            aliphatic,

        "aromaticity":
            aromaticity,

        "extinction_coefficient":
            extinction,

        "secondary_structure":
            secondary_structure,

        "flexibility":
            flexibility,

        "signal_peptide":
            signal_peptide,

        "translation":
            translation
    }


def create_protein_preview(
    protein,
    max_length=60
):

    if len(protein) <= max_length:
        return protein

    half = max_length // 2

    return (
        protein[:half]
        + "..."
        + protein[-half:]
    )
