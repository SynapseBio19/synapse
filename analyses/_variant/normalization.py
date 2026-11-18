"""Allele trimming, left normalization, contig matching and reference checks."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from analyses._variant.models import classify_variant_type, Variant, VariantAnalysisError
from analyses._variant.vcf import validate_allele


# ============================================================================
# Normalization
# ============================================================================

def trim_common_bases(
    position: int,
    reference: str,
    alternative: str,
) -> Tuple[int, str, str]:
    """
    Trim common suffix/prefix bases while preserving valid VCF-style alleles.

    The function does not left-shift through repeats; use
    left_normalize_variant() for that.
    """
    pos = position
    ref = reference.upper()
    alt = alternative.upper()

    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    return pos, ref, alt


def _reference_allele_matches(allele: str, observed: str) -> bool:
    """REF allele vs. reference bases; N in either matches any base."""
    return all(a == b or a == "N" or b == "N" for a, b in zip(allele, observed))


def left_normalize_variant(
    variant: Variant,
    reference_sequence: Optional[str] = None,
) -> Variant:
    """
    Left-normalize a simple variant.

    A reference sequence is required for repeat-aware left shifting.
    Without it, common-base trimming is still performed.
    """
    new = Variant(
        chromosome=variant.chromosome,
        position=variant.position,
        reference=variant.reference,
        alternative=variant.alternative,
        variant_id=variant.variant_id,
        quality=variant.quality,
        filter_status=variant.filter_status,
        info=dict(variant.info),
        format_keys=list(variant.format_keys),
        sample_data={
            sample: dict(values)
            for sample, values in variant.sample_data.items()
        },
        variant_type=variant.variant_type,
        gene=variant.gene,
        transcript=variant.transcript,
        region=variant.region,
        consequence=variant.consequence,
        protein_change=variant.protein_change,
        coding_change=variant.coding_change,
        allele_frequency=variant.allele_frequency,
        allele_count=variant.allele_count,
        allele_number=variant.allele_number,
        normalized=variant.normalized,
        normalization_notes=list(variant.normalization_notes),
        annotation=dict(variant.annotation),
    )

    if new.variant_type == "SYMBOLIC" or new.alternative == "*" or not validate_allele(
        new.reference, allow_symbolic=False
    ) or not validate_allele(new.alternative, allow_symbolic=False):
        new.normalization_notes.append(
            "Symbolic, spanning-deletion or non-ACGTN allele: not normalized."
        )
        return new

    pos, ref, alt = new.position, new.reference, new.alternative

    reference_checked = False
    if reference_sequence is not None:
        reference_sequence = str(reference_sequence).upper()
        observed = reference_sequence[pos - 1:pos - 1 + len(ref)]
        if len(observed) != len(ref) or not _reference_allele_matches(ref, observed):
            # Copy-on-write: the input variant's flag list is not modified (audit S16).
            _add_validation_flag(new, "REF_MISMATCH")
            new.normalization_notes.append(
                f"REF allele '{ref}' does not match the reference sequence "
                f"('{observed or 'out of range'}') at position {pos}; "
                "repeat-aware shifting was not applied."
            )
            reference_sequence = None
        else:
            reference_checked = True

    if reference_sequence is not None and len(ref) != len(alt):
        # Tan, Abecasis & Kang (2015), Bioinformatics 31:2202, Algorithm 1:
        # repeat: (1) drop an identical last base from both alleles;
        #         (2) if an allele became empty, prepend the preceding
        #             reference base to both and move one position left;
        # until neither applies. Then trim identical leading bases while
        # both alleles are longer than one base.
        changed = True
        while changed:
            changed = False
            if ref and alt and ref[-1] == alt[-1]:
                ref, alt = ref[:-1], alt[:-1]
                changed = True
            if (not ref or not alt) and pos > 1:
                previous = reference_sequence[pos - 2]
                ref, alt = previous + ref, previous + alt
                pos -= 1
                changed = True
        if (not ref or not alt) and pos == 1:
            # VCF 4.3: an indel at position 1 is anchored on the base AFTER
            # the event.
            following = reference_sequence[len(ref):len(ref) + 1]
            ref, alt = ref + following, alt + following
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref, alt = ref[1:], alt[1:]
            pos += 1
    else:
        pos, ref, alt = trim_common_bases(pos, ref, alt)

    new.position = pos
    new.reference = ref
    new.alternative = alt
    new.variant_type = classify_variant_type(ref, alt)
    new.end_position = pos + max(0, len(ref) - 1)
    new.normalized = True
    new.normalization_notes.append("Common bases trimmed.")
    if reference_checked:
        new.normalization_notes.append(
            "REF allele verified against the reference; left-aligned and "
            "parsimonious (Tan et al. 2015)."
        )
    return new


# Contig matching (audit D12). CHROM names are matched to reference names
# exactly by default. "chr_prefix" also pairs a nuclear contig with the name
# that differs only by a leading "chr" (UCSC "chr1" vs Ensembl/NCBI "1"); use
# it only when both describe the same assembly. Mitochondrial names are never
# paired, because hg19 chrM (16,571 bp) and the rCRS MT of GRCh37/GRCh38
# (16,569 bp) are different sequences. Case variants ("Chr1") are never
# folded, and CHROM itself is never renamed.
CONTIG_MATCHING_MODES = ("exact", "chr_prefix")
MITOCHONDRIAL_CONTIG_NAMES = frozenset({"M", "MT", "chrM", "chrMT"})


@dataclass
class ContigResolution:
    """How the CHROM names of a variant set map onto reference sequence names."""

    mode: str
    mapping: Dict[str, str] = field(default_factory=dict)
    aliased: Dict[str, str] = field(default_factory=dict)
    not_found: List[str] = field(default_factory=list)
    suggestions: Dict[str, str] = field(default_factory=dict)


def validate_contig_matching(contig_matching: str) -> str:
    if contig_matching not in CONTIG_MATCHING_MODES:
        raise VariantAnalysisError(
            f"Unknown contig matching mode '{contig_matching}'. "
            f"Use one of: {', '.join(CONTIG_MATCHING_MODES)}."
        )
    return contig_matching


def validate_reference_sequence(sequence: str, contig: str = "reference") -> str:
    """Raise unless *sequence* is written in the IUPAC DNA alphabet.

    Audit S22: the supplied reference was never checked, so a sequence in
    another alphabet - a peptide, a gapped alignment row, a transcript
    written with U - was compared against REF alleles base by base and
    produced mismatches that looked like data problems. Ambiguity codes are
    accepted: real reference assemblies contain N and occasionally other
    IUPAC codes.
    """
    from core.alphabets import DNA_IUPAC, normalize_sequence, symbols_outside

    invalid = symbols_outside(normalize_sequence(sequence), DNA_IUPAC)
    if invalid:
        shown = ", ".join(f"{symbol!r} at {position}" for position, symbol in invalid[:5])
        raise VariantAnalysisError(
            f"Reference sequence for '{contig}' is not DNA; found {shown}"
            + (" ..." if len(invalid) > 5 else "") + "."
        )
    return sequence


def validate_contig_references(contig_references: Mapping[str, str]) -> None:
    """Apply validate_reference_sequence() to every supplied reference."""
    for contig, sequence in contig_references.items():
        validate_reference_sequence(sequence, contig)


def _chr_prefix_partner(name: str) -> Optional[str]:
    """The name differing only by a leading "chr", or None for mitochondria."""
    if name in MITOCHONDRIAL_CONTIG_NAMES:
        return None
    partner = name[3:] if name.startswith("chr") else "chr" + name
    if not partner or partner in MITOCHONDRIAL_CONTIG_NAMES:
        return None
    return partner


def resolve_contig_names(
    contigs: Iterable[str],
    reference_names: Iterable[str],
    *,
    contig_matching: str = "exact",
) -> ContigResolution:
    """Map CHROM names onto reference names (see CONTIG_MATCHING_MODES).

    An exact name always wins. In "exact" mode a name that would match through
    the "chr" prefix is reported in ``suggestions`` and stays unmatched.
    """
    validate_contig_matching(contig_matching)
    available = set(reference_names)
    resolution = ContigResolution(mode=contig_matching)
    for contig in sorted(set(contigs)):
        if contig in available:
            resolution.mapping[contig] = contig
            continue
        partner = _chr_prefix_partner(contig)
        if partner is not None and partner in available:
            if contig_matching == "chr_prefix":
                resolution.mapping[contig] = partner
                resolution.aliased[contig] = partner
                continue
            resolution.suggestions[contig] = partner
        resolution.not_found.append(contig)
    return resolution


def _add_validation_flag(variant: Variant, flag: str) -> None:
    # Copy-on-write: normalized copies share the annotation lists of their input.
    flags = list(variant.annotation.get("validation_flags", []))
    if flag not in flags:
        flags.append(flag)
    variant.annotation["validation_flags"] = flags


def _annotate_contig_resolution(
    variant: Variant,
    resolution: ContigResolution,
    *,
    notes: bool,
) -> None:
    contig = variant.chromosome
    if contig in resolution.aliased:
        reference_name = resolution.aliased[contig]
        variant.annotation["reference_contig"] = reference_name
        _add_validation_flag(variant, "CONTIG_NAME_ALIASED")
        note = (
            f"Contig '{contig}' was matched to reference contig '{reference_name}' "
            "by the 'chr' prefix rule; CHROM is unchanged."
        )
    elif contig in resolution.not_found:
        _add_validation_flag(variant, "REFERENCE_CONTIG_NOT_FOUND")
        note = (
            f"No reference sequence for contig '{contig}': REF not checked and "
            "not left-aligned (only common bases trimmed)."
        )
        suggestion = resolution.suggestions.get(contig)
        if suggestion:
            note += (
                f" Reference contig '{suggestion}' differs only by the 'chr' prefix; "
                "use contig_matching='chr_prefix' if both describe the same assembly."
            )
    else:
        return
    if notes and note not in variant.normalization_notes:
        variant.normalization_notes.append(note)


def normalize_variants(
    variants: Sequence[Variant],
    *,
    reference_sequence: Optional[str] = None,
    contig_references: Optional[Mapping[str, str]] = None,
    contig_matching: str = "exact",
) -> List[Variant]:
    """Normalize variants.

    ``contig_references`` maps reference contig names to sequences; CHROM
    names are resolved with resolve_contig_names(). A variant whose contig has
    no reference is only trimmed (no repeat-aware shift, no REF check) and is
    flagged REFERENCE_CONTIG_NOT_FOUND; one matched through the "chr" prefix
    is flagged CONTIG_NAME_ALIASED. ``reference_sequence`` applies one
    reference to every variant (legacy use). The input variants are not
    modified.
    """
    validate_contig_matching(contig_matching)
    if reference_sequence is not None:
        validate_reference_sequence(reference_sequence)
    resolution = None
    if contig_references is not None:
        validate_contig_references(contig_references)
        resolution = resolve_contig_names(
            (v.chromosome for v in variants),
            contig_references,
            contig_matching=contig_matching,
        )

    normalized = []
    for v in variants:
        reference = reference_sequence
        if resolution is not None and contig_references is not None:
            reference_name = resolution.mapping.get(v.chromosome)
            reference = contig_references[reference_name] if reference_name else None
        result = left_normalize_variant(v, reference_sequence=reference)
        if resolution is not None:
            _annotate_contig_resolution(result, resolution, notes=True)
        normalized.append(result)
    return normalized


def check_reference_alleles(
    variants: Sequence[Variant],
    contig_references: Mapping[str, str],
    *,
    contig_matching: str = "exact",
) -> Dict[str, Any]:
    """Compare REF alleles with reference sequences (1-based VCF positions).

    CHROM names are resolved with resolve_contig_names(). Flags written to
    ``variant.annotation["validation_flags"]``: REF_MISMATCH,
    REFERENCE_CONTIG_NOT_FOUND and CONTIG_NAME_ALIASED. "not_checked" counts
    variants without a reference and variants with a symbolic or non-ACGTN
    REF; both are also reported separately.
    """
    validate_contig_references(contig_references)
    resolution = resolve_contig_names(
        (v.chromosome for v in variants),
        contig_references,
        contig_matching=contig_matching,
    )
    checked = mismatched = no_reference = not_literal = 0
    mismatches: List[Tuple[str, int, str, str]] = []
    for variant in variants:
        _annotate_contig_resolution(variant, resolution, notes=False)
        reference_name = resolution.mapping.get(variant.chromosome)
        if reference_name is None:
            no_reference += 1
            continue
        if not validate_allele(variant.reference, allow_symbolic=False):
            not_literal += 1
            continue
        reference = contig_references[reference_name].upper()
        observed = reference[variant.position - 1:variant.position - 1 + len(variant.reference)]
        checked += 1
        if len(observed) != len(variant.reference) or not _reference_allele_matches(variant.reference, observed):
            mismatched += 1
            mismatches.append(variant.key)
            _add_validation_flag(variant, "REF_MISMATCH")
    return {
        "checked": checked,
        "mismatched": mismatched,
        "not_checked": no_reference + not_literal,
        "not_checked_no_reference": no_reference,
        "not_checked_symbolic_or_non_acgtn": not_literal,
        "mismatches": mismatches,
        "contig_matching": resolution.mode,
        "contigs_not_found": list(resolution.not_found),
        "aliased_contigs": dict(resolution.aliased),
        "suggested_aliases": dict(resolution.suggestions),
    }
