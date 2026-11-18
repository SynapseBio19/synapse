"""Region annotation, sequence quality, mutation context, spectrum and hotspots."""

from __future__ import annotations

import math
from collections import Counter
from typing import List, Mapping, Sequence, Tuple

from analyses._mutation.sequence import IUPAC_COMPLEMENT, normalize_iupac_sequence
from analyses._mutation.models import (
    Mutation,
    MutationAnalysisError,
    MutationContext,
    MutationHotspot,
    MutationSpectrum,
    SequenceQuality,
)


# ============================================================================
# REGIONS
# ============================================================================

def annotate_regions(
    mutations: Sequence[Mutation],
    regions: Mapping[str, Tuple[int, int]],
) -> List[Mutation]:

    for mutation in mutations:

        for name, (
            start,
            end,
        ) in regions.items():

            if (
                mutation.position >= start
                and mutation.position <= end
            ):

                mutation.region = name
                break

    return list(mutations)


def annotate_sequence_regions(
    mutations: Sequence[Mutation],
    regions: Mapping[str, Tuple[int, int]],
) -> List[Mutation]:

    return annotate_regions(
        mutations,
        regions,
    )


# ============================================================================
# QUALITY
# ============================================================================

def sequence_complexity(
    sequence: str,
) -> float:

    dna = normalize_iupac_sequence(
        sequence
    )

    canonical = [
        base
        for base in dna
        if base in "ACGT"
    ]

    if not canonical:
        return 0.0

    counts = Counter(
        canonical
    )

    total = len(canonical)

    entropy = 0.0

    for count in counts.values():

        p = count / total

        entropy -= (
            p * math.log2(p)
        )

    # Maximum Shannon entropy for four bases = 2 bits.
    return (
        entropy / 2.0
        * 100.0
    )


def sequence_quality(
    sequence: str,
    *,
    low_complexity_threshold: float = 75.0,
) -> SequenceQuality:

    dna = normalize_iupac_sequence(
        sequence
    )

    counts = Counter(dna)

    gc = (
        100.0
        * (
            counts["G"]
            + counts["C"]
        )
        / len(dna)
    )

    n_count = counts["N"]

    ambiguous_count = sum(
        counts[base]
        for base in dna
        if base not in "ACGT"
    )

    longest = 1
    longest_base = dna[0]

    current = 1

    for i in range(
        1,
        len(dna),
    ):

        if dna[i] == dna[i - 1]:
            current += 1
        else:
            current = 1

        if current > longest:

            longest = current
            longest_base = dna[i]

    complexity = sequence_complexity(
        dna
    )

    warnings: List[str] = []

    if ambiguous_count:
        warnings.append(
            f"{ambiguous_count} ambiguous/IUPAC base(s) detected."
        )

    if longest >= 6:
        warnings.append(
            f"Long homopolymer detected: "
            f"{longest}{longest_base}."
        )

    if complexity < low_complexity_threshold:
        warnings.append(
            "Sequence shows reduced nucleotide complexity."
        )

    return SequenceQuality(
        length=len(dna),
        gc_percent=gc,
        n_count=n_count,
        ambiguous_count=ambiguous_count,
        valid=True,
        longest_homopolymer=longest,
        longest_homopolymer_base=longest_base,
        low_complexity=(
            complexity
            < low_complexity_threshold
        ),
        complexity_percent=complexity,
        warnings=warnings,
    )


# ============================================================================
# CONTEXT / SPECTRUM
# ============================================================================

def extract_mutation_context(
    mutation: Mutation,
    reference: str,
    variant: str,
    *,
    kmer_size: int = 5,
) -> MutationContext:

    if kmer_size < 1 or kmer_size % 2 == 0:
        raise MutationAnalysisError(
            "kmer_size must be a positive odd integer."
        )

    ref = normalize_iupac_sequence(
        reference
    )

    var = normalize_iupac_sequence(
        variant
    )

    radius = kmer_size // 2

    center = (
        mutation.position - 1
    )

    start = max(
        0,
        center - radius,
    )

    end = min(
        len(ref),
        center + radius + 1,
    )

    ref_context = ref[
        start:end
    ]

    if mutation.mutation_type == "Insertion":

        variant_center = min(
            len(var),
            mutation.position - 1,
        )

    else:

        variant_center = min(
            len(var) - 1,
            center,
        )

    variant_start = max(
        0,
        variant_center - radius,
    )

    variant_end = min(
        len(var),
        variant_center + radius + 1,
    )

    var_context = var[
        variant_start:variant_end
    ]

    upstream = ref[
        max(
            0,
            center - radius,
        ):
        center
    ]

    downstream = ref[
        center + 1:
        min(
            len(ref),
            center + radius + 1,
        )
    ]

    trinucleotide = None

    if (
        mutation.mutation_type == "SNP"
        and len(ref_context) >= 3
    ):

        local = center - start

        if (
            local >= 1
            and local + 1 < len(ref_context)
        ):

            trinucleotide = (
                ref_context[
                    local - 1:
                    local + 2
                ]
            )

    return MutationContext(
        position=mutation.position,
        upstream=upstream,
        downstream=downstream,
        reference_context=ref_context,
        variant_context=var_context,
        kmer_size=kmer_size,
        trinucleotide_context=trinucleotide,
    )


def build_mutation_spectrum(
    mutations: Sequence[Mutation],
) -> MutationSpectrum:

    spectrum = MutationSpectrum()

    substitution_types = (
        "A>G", "G>A",
        "C>T", "T>C",
        "A>C", "C>A",
        "A>T", "T>A",
        "C>G", "G>C",
        "G>T", "T>G",
    )

    spectrum.substitutions_12 = {
        key: 0
        for key in substitution_types
    }

    for mutation in mutations:

        if mutation.mutation_type != "SNP":
            continue

        ref = mutation.reference
        alt = mutation.alternative

        key = f"{ref}>{alt}"

        if key in spectrum.substitutions_12:
            spectrum.substitutions_12[key] += 1

        if mutation.transition is True:
            spectrum.transitions += 1

        elif mutation.transition is False:
            spectrum.transversions += 1

    # Six-class pyrimidine-oriented spectrum.
    spectrum.substitutions_6 = {
        "C>T": 0,
        "C>A": 0,
        "C>G": 0,
        "T>C": 0,
        "T>A": 0,
        "T>G": 0,
    }

    for mutation in mutations:

        if mutation.mutation_type != "SNP":
            continue

        ref = mutation.reference
        alt = mutation.alternative

        # Orient purine substitutions through complement.
        if ref in "AG":

            ref = IUPAC_COMPLEMENT[
                ref
            ]

            alt = IUPAC_COMPLEMENT[
                alt
            ]

        key = f"{ref}>{alt}"

        if key in spectrum.substitutions_6:
            spectrum.substitutions_6[key] += 1

    return spectrum


def find_mutation_hotspots(
    mutations: Sequence[Mutation],
    sequence_length: int,
    *,
    window_size: int = 20,
    minimum_mutations: int = 2,
) -> List[MutationHotspot]:

    if sequence_length < 1:
        raise MutationAnalysisError(
            "sequence_length must be positive."
        )

    if window_size < 1:
        raise MutationAnalysisError(
            "window_size must be positive."
        )

    positions = sorted(
        mutation.position
        for mutation in mutations
    )

    hotspots: List[
        MutationHotspot
    ] = []

    if not positions:
        return hotspots

    for start in range(
        1,
        sequence_length + 1,
        max(1, window_size // 2),
    ):

        end = min(
            sequence_length,
            start + window_size - 1,
        )

        count = sum(
            start <= position <= end
            for position in positions
        )

        if count >= minimum_mutations:

            density = (
                count
                / (end - start + 1)
                * 1000.0
            )

            hotspots.append(
                MutationHotspot(
                    start=start,
                    end=end,
                    mutation_count=count,
                    density_per_kb=density,
                )
            )

    # Merge overlapping windows.
    merged: List[
        MutationHotspot
    ] = []

    for hotspot in hotspots:

        if not merged:

            merged.append(hotspot)
            continue

        previous = merged[-1]

        if hotspot.start <= previous.end:

            new_start = previous.start
            new_end = max(
                previous.end,
                hotspot.end,
            )

            new_count = sum(
                start <= position <= new_end
                for position in positions
                for start in [new_start]
            )

            merged[-1] = MutationHotspot(
                start=new_start,
                end=new_end,
                mutation_count=new_count,
                density_per_kb=(
                    new_count
                    / (
                        new_end
                        - new_start
                        + 1
                    )
                    * 1000.0
                ),
            )

        else:

            merged.append(hotspot)

    return merged
