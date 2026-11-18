"""Codon resolution with IUPAC symbols and RNA translation."""

from __future__ import annotations

from data.genetic_code import CODON_TABLE
from itertools import product

from analyses._protein.constants import RNA_IUPAC


def expand_iupac_codon(codon):

    codon = (
        codon
        .upper()
        .replace("T", "U")
    )

    if len(codon) != 3:
        return []

    possibilities = []

    for base in codon:

        if base not in RNA_IUPAC:
            return []

        possibilities.append(
            RNA_IUPAC[base]
        )

    return [
        "".join(combination)
        for combination in product(
            *possibilities
        )
    ]


def rna_codon_table(genetic_code=1):
    """Codon table (RNA alphabet, "STOP" for stops) for an NCBI genetic code."""
    if genetic_code in (None, 1):
        return CODON_TABLE
    from core.seqtools import codon_dict
    return codon_dict(int(genetic_code), rna=True, stop_symbol="STOP")


def resolve_codon(codon, codon_table=None):
    codon_table = CODON_TABLE if codon_table is None else codon_table

    codon = (
        codon
        .upper()
        .replace("T", "U")
    )

    expanded = expand_iupac_codon(
        codon
    )

    if not expanded:

        return {
            "codon": codon,
            "status": "invalid",
            "amino_acid": None,
            "possibilities": [],
            "expanded_codons": [],
            "contains_stop": False,
        }

    translations = set()

    for concrete_codon in expanded:

        amino_acid = codon_table.get(
            concrete_codon
        )

        if amino_acid is not None:

            translations.add(
                amino_acid
            )

    if not translations:

        return {
            "codon": codon,
            "status": "invalid",
            "amino_acid": None,
            "possibilities": [],
            "expanded_codons": expanded,
            "contains_stop": False,
        }

    contains_stop = (
        "STOP" in translations
    )

    amino_acids = sorted(
        amino_acid
        for amino_acid in translations
        if amino_acid != "STOP"
    )

    if len(translations) == 1:

        result = next(
            iter(translations)
        )

        if result == "STOP":

            return {
                "codon": codon,
                "status": "stop",
                "amino_acid": None,
                "possibilities": ["STOP"],
                "expanded_codons": expanded,
                "contains_stop": True,
            }

        return {
            "codon": codon,
            "status": "resolved",
            "amino_acid": result,
            "possibilities": [result],
            "expanded_codons": expanded,
            "contains_stop": False,
        }

    if not contains_stop:

        return {
            "codon": codon,
            "status":
                "ambiguous_amino_acid",
            "amino_acid": None,
            "possibilities":
                amino_acids,
            "expanded_codons": expanded,
            "contains_stop": False,
        }

    return {
        "codon": codon,
        "status": "ambiguous_stop",
        "amino_acid": None,
        "possibilities":
            amino_acids + ["STOP"],
        "expanded_codons": expanded,
        "contains_stop": True,
    }


# ============================================================
# TRANSLATION
# ============================================================


def translate_protein(rna):

    if not isinstance(rna, str):
        return ""

    rna = (
        "".join(rna.upper().split())
        .replace("T", "U")
    )

    protein = ""

    for i in range(
        0,
        len(rna) - 2,
        3
    ):

        codon = rna[
            i:i + 3
        ]

        result = resolve_codon(
            codon
        )

        if result["status"] == "stop":
            break

        if result["status"] == "resolved":

            protein += (
                result["amino_acid"]
            )

        else:

            protein += "?"

    return protein


def analyze_translation(rna, initiator_start=False, genetic_code=1):
    """Translate RNA in frame +1 with the standard genetic code.

    initiator_start=True states that the first codon is the selected start
    codon of an ORF. If it is an initiation codon of the genetic code
    (analyses.orf.is_initiation_codon; table 1: AUG, CUG, GUG, UUG) it is
    decoded as methionine, because initiator tRNA inserts (formyl-)Met at
    the start codon (audits P0-8, D1); any other first codon is read with
    the elongation code. For arbitrary sequence windows leave it False.
    genetic_code selects the NCBI translation table (default 1, standard).
    """
    from analyses.orf import is_initiation_codon

    if not isinstance(rna, str):
        return {
            "rna": "",
            "reading_frame": "+1",
            "protein": "",
            "protein_length": 0,
            "codon_count_available": 0,
            "codons_processed": 0,
            "codon_records": [],
            "start_codon": None,
            "starts_with_aug": False,
            "stop_found": False,
            "stop_codon": None,
            "stop_position": None,
            "trailing_nucleotides": 0,
            "ambiguous_codon_count": 0,
            "deterministic_ambiguous_count": 0,
            "uncertain_stop_count": 0,
            "invalid_codon_count": 0,
        }

    rna = (
        "".join(rna.upper().split())
        .replace("T", "U")
    )

    available_codons = (
        len(rna) // 3
    )

    trailing_nucleotides = (
        len(rna) % 3
    )

    codon_table = rna_codon_table(genetic_code)
    records = []

    protein = ""

    stop_found = False
    stop_codon = None
    stop_position = None

    ambiguous_codons = 0
    deterministic_ambiguous = 0
    uncertain_stop_codons = 0
    invalid_codons = 0

    for codon_index in range(
        available_codons
    ):

        start = codon_index * 3
        end = start + 3

        codon = rna[
            start:end
        ]

        result = resolve_codon(codon, codon_table)

        initiator = False
        if (
            initiator_start
            and codon_index == 0
            and is_initiation_codon(codon, genetic_code)
        ):
            result = dict(result)
            result["status"] = "resolved"
            result["amino_acid"] = "M"
            result["possibilities"] = ["M"]
            initiator = True

        iupac_ambiguous = any(
            base not in "AUCG"
            for base in codon
        )

        if iupac_ambiguous:
            ambiguous_codons += 1

        record = {
            "number":
                codon_index + 1,
            "codon":
                codon,
            "start":
                start + 1,
            "end":
                end,
            "status":
                result["status"],
            "amino_acid":
                result["amino_acid"],
            "possibilities":
                result["possibilities"],
            "contains_stop":
                result["contains_stop"],
            "iupac_ambiguous":
                iupac_ambiguous,
            "initiator":
                initiator,
        }

        records.append(record)

        if result["status"] == "stop":

            stop_found = True
            stop_codon = codon
            stop_position = start + 1
            break

        if result["status"] == "resolved":

            protein += (
                result["amino_acid"]
            )

            if iupac_ambiguous:
                deterministic_ambiguous += 1

        elif result["status"] == "ambiguous_amino_acid":

            protein += "?"

        elif result["status"] == "ambiguous_stop":

            protein += "?"
            uncertain_stop_codons += 1

        else:

            protein += "?"
            invalid_codons += 1

    starts_with_aug = (
        len(rna) >= 3
        and rna[:3] == "AUG"
    )

    return {
        "rna":
            rna,

        "reading_frame":
            "+1",

        "genetic_code":

            int(genetic_code or 1),

        "protein":
            protein,

        "protein_length":
            len(protein),

        "codon_count_available":
            available_codons,

        "codons_processed":
            len(records),

        "codon_records":
            records,

        "start_codon":
            rna[:3]
            if len(rna) >= 3
            else None,

        "starts_with_aug":
            starts_with_aug,

        "stop_found":
            stop_found,

        "stop_codon":
            stop_codon,

        "stop_position":
            stop_position,

        "trailing_nucleotides":
            trailing_nucleotides,

        "ambiguous_codon_count":
            ambiguous_codons,

        "deterministic_ambiguous_count":
            deterministic_ambiguous,

        "uncertain_stop_count":
            uncertain_stop_codons,

        "invalid_codon_count":
            invalid_codons,

        "initiator_start":
            bool(initiator_start),
    }
