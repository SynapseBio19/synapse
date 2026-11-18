"""HGVS c. and p. descriptions (3' rule, dup, delins, alleles)."""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from analyses._mutation.models import Mutation


# ============================================================================
# HGVS-LIKE NOTATION
# ============================================================================

def mutation_hgvs_c(
    mutation: Mutation,
    reference_sequence: Optional[str] = None,
) -> str:
    """HGVS DNA description of one mutation ("c." + HGVS syntax).

    Forms (HGVS recommendations, den Dunnen et al. 2016): ``c.7C>T``;
    ``c.16del`` / ``c.10_12del`` (the deleted bases are not repeated);
    ``c.10_11insTT`` between the two flanking bases; ``c.6dup`` /
    ``c.12_13dup`` when the inserted bases repeat the bases just before
    the insertion site; ``c.11_12delinsCCC`` for other replacements.

    With ``reference_sequence`` (the ungapped reference the mutation was
    called on) deletions and insertions are moved to their most 3'
    position in a repeat (the HGVS 3' rule) and duplications are
    recognized; without it the positions are written as detected. Checked
    against Mutalyzer 3 (``description_extract``) on synthetic cases, see
    docs/METHODS.md. Position 1 is the first base of the supplied sequence,
    which equals HGVS c. numbering only if that sequence is the CDS starting
    at the A of the start codon.
    """

    sequence = (
        reference_sequence.upper()
        if reference_sequence
        else None
    )
    reference = mutation.reference
    alternative = mutation.alternative

    if mutation.mutation_type == "Insertion":
        # Inserted before `position`, after `position - 1`.
        after = mutation.position
        inserted = alternative
        if sequence and inserted:
            while (
                after <= len(sequence)
                and inserted[0] == sequence[after - 1]
            ):
                inserted = inserted[1:] + inserted[0]
                after += 1
            length = len(inserted)
            if (
                after - 1 >= length
                and sequence[after - 1 - length:after - 1] == inserted
            ):
                first, last = after - length, after - 1
                return (
                    f"c.{first}dup"
                    if first == last
                    else f"c.{first}_{last}dup"
                )
        before = after - 1
        # HGVS has no position 0: the base before c.1 is c.-1.
        return f"c.{before if before >= 1 else -1}_{after}ins{inserted}"

    if mutation.mutation_type == "Deletion":
        first, last = mutation.position, mutation.end_position
        if (
            sequence
            and sequence[first - 1:last] == reference.upper()
        ):
            while (
                last < len(sequence)
                and sequence[first - 1] == sequence[last]
            ):
                first += 1
                last += 1
        return (
            f"c.{first}del"
            if first == last
            else f"c.{first}_{last}del"
        )

    if len(reference) == 1 and len(alternative) == 1:
        return f"c.{mutation.position}{reference}>{alternative}"

    first, last = mutation.position, mutation.end_position
    location = f"{first}" if first == last else f"{first}_{last}"
    return f"c.{location}delins{alternative}"


def _reference_span(mutation: Mutation) -> Tuple[int, int]:
    """1-based reference interval; an insertion before p spans (p, p - 1)."""
    if mutation.mutation_type == "Insertion":
        return mutation.position, mutation.position - 1
    return mutation.position, max(mutation.position, mutation.end_position)


def hgvs_c_allele(
    mutations: Sequence[Mutation],
    reference_sequence: str,
) -> Optional[str]:
    """HGVS description of all mutations as one allele, e.g. c.[7C>T;20_21delinsCC].

    Events that touch (no unchanged reference base between them) are one
    change in HGVS and are merged into a deletion-insertion; a single
    event uses mutation_hgvs_c() with the 3' rule. One change is written
    without brackets. None when there are no mutations.
    """
    if not mutations:
        return None
    sequence = reference_sequence.upper()
    ordered = sorted(mutations, key=lambda item: (_reference_span(item)[0], item.mutation_id))
    groups: List[List[Mutation]] = []
    for mutation in ordered:
        start, _ = _reference_span(mutation)
        if groups and start <= max(_reference_span(item)[1] for item in groups[-1]) + 1:
            groups[-1].append(mutation)
        else:
            groups.append([mutation])
    descriptions = []
    for group in groups:
        if len(group) == 1:
            descriptions.append(mutation_hgvs_c(group[0], sequence)[2:])
            continue
        first = min(_reference_span(item)[0] for item in group)
        last = max(_reference_span(item)[1] for item in group)
        observed, cursor = [], first
        for item in group:
            start, end = _reference_span(item)
            observed.append(sequence[cursor - 1:start - 1])
            observed.append(item.alternative)
            cursor = max(cursor, end + 1)
        observed.append(sequence[cursor - 1:last])
        location = f"{first}" if first == last else f"{first}_{last}"
        descriptions.append(f"{location}delins{''.join(observed)}")
    if len(descriptions) == 1:
        return f"c.{descriptions[0]}"
    return "c.[" + ";".join(descriptions) + "]"


def mutation_hgvs_p(
    mutation: Mutation,
) -> Optional[str]:
    """Predicted protein consequence in HGVS style, one-letter codes.

    The consequence is inferred from the DNA, not observed, so it is
    written in parentheses (HGVS recommendations, den Dunnen et al. 2016):
    p.(R12C) missense, p.(R12=) synonymous, p.(W24*) nonsense, p.(M1?)
    for a change of the initiator codon (the effect on initiation is
    unknown) and p.(*110Qext*?) for a lost stop codon whose new stop is not
    determined. Positions count codons of the supplied sequence, which is
    HGVS p. numbering only if that sequence is the CDS. Not checked
    against an HGVS validator. None when no codon change is known.
    """

    if (
        mutation.reference_amino_acid
        is None
        or mutation.alternative_amino_acid
        is None
        or mutation.codon_index is None
    ):
        return None

    ref = mutation.reference_amino_acid
    alt = mutation.alternative_amino_acid
    index = mutation.codon_index

    if ref == alt:
        return f"p.({ref}{index}=)"

    if index == 1 and ref == "M":
        return "p.(M1?)"

    if ref == "*":
        return f"p.(*{index}{alt}ext*?)"

    return f"p.({ref}{index}{alt})"
