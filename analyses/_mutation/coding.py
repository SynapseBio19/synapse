"""Codon and protein consequences of mutations and the coding-frame check."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from analyses._mutation.sequence import (
    CODON_TABLE,
    normalize_dna,
    normalize_iupac_sequence,
    reverse_complement_iupac,
)
from analyses._mutation.models import CodingAnalysis, Mutation, MutationAnalysisError


def classify_coding_substitution(
    reference_codon: str,
    alternative_codon: str,
    codon_index: Optional[int] = None,
) -> str:
    """Classify a codon substitution (standard genetic code).

    ``codon_index`` is the 1-based codon number within the CDS. "Start-loss"
    applies only to codon 1 (the initiator); a Met->X change at any other
    codon is an ordinary missense change (audit P0-6). When ``codon_index``
    is None the position is unknown and Met->X is reported as "Missense".
    """

    ref = normalize_dna(
        reference_codon,
        allow_n=False,
    )

    alt = normalize_dna(
        alternative_codon,
        allow_n=False,
    )

    if len(ref) != 3 or len(alt) != 3:
        raise MutationAnalysisError(
            "Codons must contain exactly 3 bases."
        )

    ref_aa = CODON_TABLE.get(
        ref,
        "X",
    )

    alt_aa = CODON_TABLE.get(
        alt,
        "X",
    )

    if ref_aa == alt_aa:
        return "Synonymous"

    if ref_aa != "*" and alt_aa == "*":
        return "Nonsense"

    if ref_aa == "*" and alt_aa != "*":
        return "Stop-loss"

    if codon_index == 1 and ref_aa == "M" and alt_aa != "M":
        return "Start-loss"

    return "Missense"


def _translate_canonical(
    sequence: str,
) -> str:

    dna = normalize_dna(
        sequence,
        allow_n=False,
    )

    protein = []

    for i in range(
        0,
        len(dna) - 2,
        3,
    ):

        protein.append(
            CODON_TABLE.get(
                dna[i:i + 3],
                "X",
            )
        )

    return "".join(protein)


def _first_difference_position(reference: str, alternative: str) -> Optional[int]:
    """Return the first one-based coordinate where two sequences differ."""
    for index, (ref_base, alt_base) in enumerate(zip(reference, alternative), 1):
        if ref_base != alt_base:
            return index
    if len(reference) != len(alternative):
        return min(len(reference), len(alternative)) + 1
    return None


def analyze_coding_consequences(
    reference_cds: str,
    alternative_cds: str,
) -> CodingAnalysis:
    """
    Analyze coding consequences without treating an indel as a cascade of
    independent downstream amino-acid substitutions.

    Rules:
    - Equal coding length: compare complete codons one-to-one.
    - Length difference divisible by 3: classify the event as an in-frame
      indel; do not count shifted downstream codons as missense mutations.
    - Length difference not divisible by 3: classify the event as a frameshift;
      downstream amino-acid differences are not counted as individual
      missense mutations.
    """
    ref = normalize_iupac_sequence(reference_cds)
    alt = normalize_iupac_sequence(alternative_cds)

    if len(ref) < 3:
        raise MutationAnalysisError(
            "Reference CDS must contain at least one complete codon."
        )

    # Ambiguous IUPAC symbols are valid nucleotide input, but they do not
    # define a unique amino acid. Preserve the sequence-level information
    # without inventing a protein consequence.
    canonical_ref = set(ref) <= set("ACGT")
    canonical_alt = set(alt) <= set("ACGT")

    # A reverse-complement sequence represents the same nucleotide sequence
    # in the opposite orientation. For coding-consequence analysis, normalize
    # that orientation before comparing codons. Otherwise a pure orientation
    # change would be incorrectly reported as dozens of missense changes.
    ref_protein = (
        _translate_canonical(ref)
        if canonical_ref
        else ""
    )
    alt_protein = (
        _translate_canonical(alt)
        if canonical_alt
        else ""
    )

    reverse_complement_match = (
        len(ref) == len(alt)
        and alt == reverse_complement_iupac(ref)
        and alt != ref
    )

    if reverse_complement_match and canonical_ref and canonical_alt:
        coding = CodingAnalysis(
            reference_protein=ref_protein,
            alternative_protein=ref_protein,
            reference_coding_length=len(ref),
            alternative_coding_length=len(alt),
            protein_length_change=0,
            notes=[
                "Alternative CDS is the reverse complement of the reference CDS.",
                "No coding consequences were reported because the difference is orientation-only.",
            ],
        )
        return coding

    if not (canonical_ref and canonical_alt):
        return CodingAnalysis(
            reference_protein="",
            alternative_protein="",
            reference_coding_length=len(ref),
            alternative_coding_length=len(alt),
            protein_length_change=0,
            notes=[
                "Coding consequence: Unknown.",
                "Exact amino-acid consequences were not inferred because "
                "ambiguous IUPAC bases are present.",
                "A codon containing an ambiguous nucleotide cannot be assigned "
                "a unique amino acid without resolving the ambiguity.",
            ],
        )

    coding = CodingAnalysis(
        reference_protein=ref_protein,
        alternative_protein=alt_protein,
        reference_coding_length=len(ref),
        alternative_coding_length=len(alt),
        protein_length_change=len(alt_protein) - len(ref_protein),
    )

    delta = len(alt) - len(ref)

    # No length change: true codon-by-codon substitution analysis.
    if delta == 0:
        complete_codons = min(len(ref), len(alt)) // 3

        for codon_zero in range(complete_codons):
            start = codon_zero * 3
            ref_codon = ref[start:start + 3]
            alt_codon = alt[start:start + 3]

            if ref_codon == alt_codon:
                continue

            effect = classify_coding_substitution(
                ref_codon,
                alt_codon,
                codon_index=codon_zero + 1,
            )

            coding.affected_codons.append(codon_zero + 1)

            if effect == "Synonymous":
                coding.synonymous_count += 1
            elif effect == "Missense":
                coding.missense_count += 1
            elif effect == "Nonsense":
                coding.nonsense_count += 1
                if coding.premature_stop_position is None:
                    coding.premature_stop_position = codon_zero + 1
            elif effect == "Start-loss":
                coding.start_loss_count += 1
            elif effect == "Stop-loss":
                coding.stop_loss_count += 1

        return coding

    # Any non-zero length change is structurally an indel event here.
    first_diff = _first_difference_position(ref, alt)
    if first_diff is None:
        first_diff = min(len(ref), len(alt)) + 1

    codon_index = (first_diff - 1) // 3 + 1
    coding.affected_codons.append(codon_index)

    if delta % 3 == 0:
        coding.notes.append(
            f"In-frame indel: coding length changed by {delta:+d} bp."
        )
        return coding

    coding.frameshift = True
    coding.frameshift_position = first_diff
    coding.notes.append(
        f"Frameshift: coding length changed by {delta:+d} bp."
    )

    # A stop after a frameshift is descriptive only; it is not counted as a
    # separate nonsense substitution.
    frame_aa_index = max(0, (first_diff - 1) // 3)
    for aa_index in range(frame_aa_index, len(alt_protein)):
        if alt_protein[aa_index] == "*":
            coding.premature_stop_position = aa_index + 1
            break

    return coding


def annotate_mutations_with_coding_consequences(
    mutations: List[Mutation],
    reference_cds: str,
    alternative_cds: str,
) -> CodingAnalysis:
    """
    Add coding annotations to detected events.

    Important:
    - SNPs/substitutions are classified codon-by-codon only when both
      sequences have identical length.
    - Insertions/deletions are classified structurally as either in-frame
      or frameshift events.
    - IUPAC ambiguity is accepted at sequence level, but an exact amino-acid
      consequence is not invented for ambiguous codons.
    """
    ref = normalize_iupac_sequence(reference_cds)
    alt = normalize_iupac_sequence(alternative_cds)

    canonical_ref = set(ref) <= set("ACGT")
    canonical_alt = set(alt) <= set("ACGT")

    # Structural annotations are always safe and do not require canonical CDS.
    for mutation in mutations:
        mutation.effect = None
        mutation.consequence = None
        mutation.frameshift = False
        mutation.premature_stop = False

        mutation.codon_index = (
            (max(1, mutation.position) - 1) // 3 + 1
        )
        mutation.codon_position = (
            (max(1, mutation.position) - 1) % 3 + 1
        )

        if mutation.mutation_type in {"Insertion", "Deletion"}:
            if mutation.length_change % 3 == 0:
                mutation.effect = "In-frame indel"
                mutation.consequence = "In-frame indel"
            else:
                mutation.effect = "Frameshift"
                mutation.consequence = "Frameshift"
                mutation.frameshift = True
            continue

        if mutation.mutation_type == "Complex":
            if mutation.length_change % 3 == 0:
                mutation.effect = "In-frame indel"
                mutation.consequence = "In-frame indel"
            else:
                mutation.effect = "Frameshift"
                mutation.consequence = "Frameshift"
                mutation.frameshift = True

    # Ambiguous bases: preserve sequence-level detection but do not fabricate
    # exact protein consequences.
    if not (canonical_ref and canonical_alt):
        coding = CodingAnalysis(
            reference_protein="",
            alternative_protein="",
            reference_coding_length=len(ref),
            alternative_coding_length=len(alt),
            protein_length_change=0,
            notes=[
                "Coding consequence: Unknown.",
                "Exact amino-acid consequences were not inferred because "
                "ambiguous IUPAC bases are present.",
                "A codon containing an ambiguous nucleotide cannot be assigned "
                "a unique amino acid without resolving the ambiguity.",
            ],
        )
        for mutation in mutations:
            if mutation.mutation_type in {"SNP", "Substitution"}:
                mutation.effect = "Unknown"
                mutation.consequence = "Unknown"
                mutation.notes.append(
                    "Ambiguous IUPAC base prevents exact codon consequence."
                )
        return coding

    coding = analyze_coding_consequences(
        ref,
        alt,
    )

    # Only equal-length sequences support direct SNP codon comparison.
    # For indel-containing comparisons, this deliberately avoids the old
    # downstream-missense cascade.
    if len(ref) != len(alt):
        return coding

    complete_codons = len(ref) // 3

    for mutation in mutations:
        if mutation.mutation_type not in {"SNP", "Substitution"}:
            continue

        codon_zero = (mutation.position - 1) // 3
        codon_start = codon_zero * 3

        if codon_zero >= complete_codons:
            mutation.effect = "Unknown"
            mutation.consequence = "Unknown"
            mutation.notes.append(
                "Mutation lies in an incomplete terminal codon."
            )
            continue

        ref_codon = ref[codon_start:codon_start + 3]
        alt_codon = alt[codon_start:codon_start + 3]

        if len(ref_codon) != 3 or len(alt_codon) != 3:
            mutation.effect = "Unknown"
            mutation.consequence = "Unknown"
            continue

        mutation.codon_index = codon_zero + 1
        mutation.codon_position = (mutation.position - 1) % 3 + 1
        mutation.reference_codon = ref_codon
        mutation.alternative_codon = alt_codon
        mutation.reference_amino_acid = CODON_TABLE[ref_codon]
        mutation.alternative_amino_acid = CODON_TABLE[alt_codon]

        effect = classify_coding_substitution(
            ref_codon,
            alt_codon,
            codon_index=codon_zero + 1,
        )
        mutation.effect = effect
        mutation.consequence = effect

        if effect == "Nonsense":
            mutation.premature_stop = True

    return coding

CODING_FRAME_ASSUMPTION = (
    "Coding consequences assume that position 1 of the reference is the first "
    "base of codon 1 (CDS coordinates). For a genomic or non-coding sequence "
    "this assumption is wrong and amino-acid effects are meaningless."
)


def assess_coding_frame(cds: str) -> List[str]:
    """Return warnings when a sequence does not look like a CDS in frame +1.

    These checks are plausibility checks, not CDS validation: a sequence can
    pass all of them and still not be the real coding sequence.
    """
    seq = normalize_iupac_sequence(cds)
    warnings: List[str] = []
    if len(seq) % 3:
        warnings.append(
            f"Length {len(seq)} is not a multiple of 3; the last "
            f"{len(seq) % 3} base(s) do not form a codon."
        )
    if len(seq) >= 3 and seq[:3] not in {"ATG", "GTG", "TTG", "CTG"}:
        warnings.append(
            f"Does not begin with a start codon ({seq[:3]}); "
            "frame +1 may not be the coding frame."
        )
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    internal_stops = [
        index + 1
        for index, codon in enumerate(codons[:-1])
        if codon in {"TAA", "TAG", "TGA"}
    ]
    if internal_stops:
        shown = ", ".join(str(i) for i in internal_stops[:5])
        warnings.append(
            f"{len(internal_stops)} in-frame stop codon(s) before the end "
            f"(codon {shown}{', ...' if len(internal_stops) > 5 else ''}); "
            "reading frame +1 is unlikely to be the coding frame."
        )
    if codons and codons[-1] not in {"TAA", "TAG", "TGA"}:
        warnings.append("Does not end with a stop codon (possibly partial CDS).")
    return warnings


def orf_reference_cds(orf: Any, reference: str, variant: str) -> Tuple[str, Optional[str], str]:
    """Derive reference/variant CDS from a pipeline ORF.

    Returns (reference_cds, variant_cds or None, description). The variant
    CDS is only derived when the variant has the reference length (pure
    substitutions), because only then do ORF coordinates carry over directly.
    """
    from core.orf_context import resolve_orf_region

    region = resolve_orf_region(orf, reference)
    # An ORF without a sequence gives no region (None); before 2026-11-12
    # that raised AttributeError instead of this error (D22).
    if (region is None or not region.validated or region.start_1based is None
            or region.end_1based is None):
        raise MutationAnalysisError(
            "The selected ORF could not be located on the reference sequence."
        )
    ref_cds = region.coding_sequence
    description = (
        f"{region.name}: {region.start_1based}-{region.end_1based} "
        f"({region.strand} strand)"
    )
    var_cds = None
    if len(variant) == len(reference):
        forward = variant[region.start_1based - 1:region.end_1based]
        var_cds = reverse_complement_iupac(forward) if region.strand == "-" else forward
    return ref_cds, var_cds, description


def compare_proteins(
    reference_protein: str,
    alternative_protein: str,
) -> Dict[str, Any]:

    ref = str(
        reference_protein
    )

    alt = str(
        alternative_protein
    )

    length = min(
        len(ref),
        len(alt),
    )

    mismatches = []

    for i in range(length):

        if ref[i] != alt[i]:

            mismatches.append(
                {
                    "position": i + 1,
                    "reference": ref[i],
                    "alternative": alt[i],
                }
            )

    return {
        "reference_length": len(ref),
        "alternative_length": len(alt),
        "length_change": len(alt) - len(ref),
        "identity_percent": (
            100.0
            if max(len(ref), len(alt)) == 0
            else (
                100.0
                * (
                    len(ref)
                    - len(mismatches)
                )
                / max(
                    len(ref),
                    len(alt),
                )
            )
        ),
        "mismatches": mismatches,
    }
