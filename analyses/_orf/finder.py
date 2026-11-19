"""The ORF search: complete, 5'-partial and origin-spanning ORFs."""

from __future__ import annotations

from analyses.sequence import get_reverse_complement

from analyses._orf.codes import (
    analyze_codon_identity,
    codon_table_for,
    DEFAULT_MIN_ORF_LENGTH,
    DEFAULT_RANKING,
    DEFAULT_REPORT_NESTED,
    ORF_RANKING_MODES,
    resolve_start_codons,
    start_codons_for,
    _stop_codons_of,
)
from analyses._orf.features import (
    analyze_codon_usage,
    analyze_orf_ambiguity,
    analyze_start_context,
    calculate_cds_quality,
    calculate_codon_diversity,
    calculate_orf_heuristic_score,
)


def _orf_sort_key(orf, ranking):
    if ranking == "heuristic":
        return (
            orf["heuristic_score"]["score"],
            orf["cds_quality"]["score"],
            orf["length"],
            orf["status"] == "Complete",
            orf["start_type"] == "Canonical",
        )
    # "length": longest ORF first (ORFfinder-style); ties broken by
    # completeness, canonical start and then the heuristic score.
    return (
        orf["length"],
        orf["status"] == "Complete",
        orf["start_type"] == "Canonical",
        orf["heuristic_score"]["score"],
    )


def _keep_longest_per_stop(orfs):
    """Keep only the outermost ORF for each (strand, 3' end, status).

    Strand and 3' end fix the reading frame on a linear sequence, so the frame
    adds nothing there; on a circle whose length is not a multiple of 3 the
    frame changes across the origin, and comparing frames would keep two ORFs
    that end at the same stop.
    """
    best = {}
    for orf in orfs:
        # An ORF across the origin ends at stop_on_circle; it competes with
        # the ORFs ending at that base.
        three_prime = (orf.get("stop_on_circle", orf["stop"])
                       if orf["strand"] == "+" else orf["start"])
        key = (orf["strand"], three_prime, orf["status"])
        current = best.get(key)
        if current is None or orf["length"] > current["length"]:
            best[key] = orf
    return list(best.values())


# Statuses of ORFs that do not begin at a start codon inside the sequence.
MISSING_START_STATUSES = ("Partial - Missing Start", "Partial - Missing Start and Stop")


def orf_has_start_codon(orf_or_status):
    """False for a 5'-partial ORF, whose first codon is not a start codon."""
    status = orf_or_status.get("status") if isinstance(orf_or_status, dict) else orf_or_status
    return status not in MISSING_START_STATUSES


def _five_prime_partial_orfs(dna, *, min_length, start_codons, stop_codons, codon_table,
                             genetic_code, include_partial):
    """Open frames that run in from the 5' end of a linear sequence (P1-9).

    In each frame the segment from the first codon to the first stop is
    open but has no start codon in the sequence: the start lies upstream,
    outside it. It is reported as "Partial - Missing Start" when it reaches
    min_length; with include_partial also a frame without any stop, as
    "Partial - Missing Start and Stop". A frame whose first codon is itself
    a start codon gives an ordinary complete ORF and no partial one. Complete
    ORFs inside the segment are kept: the per-stop rule compares ORFs of the
    same status only.
    """
    records = []
    length = len(dna)
    for strand, sequence in (("+", dna), ("-", get_reverse_complement(dna))):
        for offset in range(3):
            if offset + 3 > len(sequence):
                continue
            first_codon = sequence[offset:offset + 3]
            if analyze_codon_identity(first_codon, start_codons) == "definite":
                continue
            frame_number = offset + 1 if strand == "+" else -(offset + 1)
            possible_stop_sites = []
            stop_codon = None
            position = offset
            while position <= len(sequence) - 3:
                codon = sequence[position:position + 3]
                status = analyze_codon_identity(codon, stop_codons)
                if status == "definite":
                    stop_codon = codon
                    break
                if status == "possible":
                    possible_stop_sites.append({"codon": codon, "position": position + 1})
                position += 3
            if stop_codon is not None:
                end = position + 3
                status_text = "Partial - Missing Start"
            elif include_partial:
                end = offset + ((len(sequence) - offset) // 3) * 3
                status_text = "Partial - Missing Start and Stop"
            else:
                continue
            orf_sequence = sequence[offset:end]
            if len(orf_sequence) < max(3, min_length):
                continue
            if strand == "+":
                genomic_start, genomic_stop = offset, end - 1
            else:
                genomic_start, genomic_stop = length - end, length - offset - 1
            ambiguity = analyze_orf_ambiguity(orf_sequence, codon_table=codon_table)
            codon_count = len(orf_sequence) // 3
            orf = {
                "genetic_code": genetic_code,
                "strand": strand,
                "frame": frame_number,
                "start": genomic_start,
                "stop": genomic_stop,
                "length": len(orf_sequence),
                "codon_count": codon_count,
                "protein_length": codon_count - 1 if stop_codon else codon_count,
                "start_codon": first_codon,
                "start_type": "None (5' partial)",
                "start_status": "absent",
                "stop_codon": stop_codon,
                "status": status_text,
                "sequence": orf_sequence,
                "possible_stop_sites": possible_stop_sites,
                "ambiguous_codon_count": ambiguity["ambiguous_codon_count"],
                "ambiguous_codons": ambiguity["ambiguous_codons"],
                "deterministic_ambiguous_count": ambiguity["deterministic_ambiguous_count"],
                "uncertain_codon_count": ambiguity["uncertain_codon_count"],
                "possible_stop_codon_count": ambiguity["possible_stop_codon_count"],
                "translation_safe": ambiguity["translation_safe"],
                "start_context": analyze_start_context(sequence, offset),
                "codon_usage": analyze_codon_usage(orf_sequence, codon_table=codon_table),
                "codon_diversity": calculate_codon_diversity(orf_sequence, codon_table=codon_table),
            }
            orf["cds_quality"] = calculate_cds_quality(orf)
            orf["heuristic_score"] = calculate_orf_heuristic_score(orf)
            records.append(orf)
    return records


def _origin_spanning_orfs(dna, **search):
    """Complete ORFs that run across the origin of a circular sequence (P1-9).

    The sequence is searched once more as dna + dna[:n-1], which lets a frame
    continue across the origin. An ORF from that search is kept when it
    starts in the first copy and ends in the second, and it is not longer
    than the circle. Coordinates stay unrolled, so start <= stop and the
    length, sequence and interval arithmetic hold; "stop_on_circle" gives the
    last base on the circle and "wraps_origin" marks the ORF.
    """
    length = len(dna)
    if length < 6:
        return []
    doubled = dna + dna[:length - 1]
    wrapped = []
    for orf in find_orfs(doubled, include_partial=False, nested=True, **search):
        if orf["start"] < length <= orf["stop"] and orf["length"] <= length:
            orf["wraps_origin"] = True
            orf["stop_on_circle"] = orf["stop"] - length
            # the frame of its start codon on the circle (strand-relative)
            orf["frame"] = (orf["start"] % 3 + 1 if orf["strand"] == "+"
                            else -((length - 1 - orf["stop_on_circle"]) % 3 + 1))
            wrapped.append(orf)
    return wrapped


def orf_end_on_sequence(orf):
    """1-based last base of an ORF on its sequence (across the origin: on the circle)."""
    return orf.get("stop_on_circle", orf["stop"]) + 1


def find_orfs(
    dna,
    include_partial=True,
    min_length=DEFAULT_MIN_ORF_LENGTH,
    start_codons=None,
    stop_codons=None,
    start_codon_policy=None,
    nested=DEFAULT_REPORT_NESTED,
    ranking=DEFAULT_RANKING,
    genetic_code=1,
    topology="linear",
    include_5prime_partial=False,
):
    """Find open reading frames on both strands.

    Parameters
    ----------
    min_length:
        Minimum ORF length in nucleotides, stop codon included.
    start_codon_policy:
        "atg_only" (default) or "atg_and_alternative" (ATG/GTG/TTG). Ignored
        when ``start_codons`` is given explicitly.
    nested:
        False (default): report only the longest ORF per stop codon.
        True: report every in-frame start as its own (nested) ORF.
    ranking:
        "length" (default, longest first) or "heuristic" (legacy
        heuristic score).
    topology:
        "circular" adds the ORFs that run across the origin
        (wraps_origin, stop_on_circle) and drops the partial ORFs, which on
        a circle are artifacts of the arbitrary origin.
    include_5prime_partial:
        Also report open frames entering from the 5' end without a start
        codon ("Partial - Missing Start"; linear sequences only).

    Coordinates in the result are 0-based and inclusive on the forward
    strand ("start" <= "stop"), for both strands.

    The sequence is uppercased with whitespace removed (soft-masked input
    is searched like any other); symbols outside the IUPAC DNA alphabet,
    RNA's U included, raise ValueError (audit D14-F8). ``core.api``
    converts RNA explicitly before calling this.
    """

    from core.alphabets import DNA_IUPAC, normalize_sequence, symbols_outside

    dna = normalize_sequence(dna)
    invalid = symbols_outside(dna, DNA_IUPAC)
    if invalid:
        shown = ", ".join(f"{symbol!r} at {position}" for position, symbol in invalid[:5])
        raise ValueError(
            f"ORF search needs a DNA sequence; found {shown}"
            + (" ..." if len(invalid) > 5 else "") + "."
        )

    if ranking not in ORF_RANKING_MODES:
        raise ValueError(
            f"ranking must be one of {ORF_RANKING_MODES}."
        )

    codon_table = codon_table_for(genetic_code)
    if start_codons is None and (genetic_code not in (None, 1) or start_codon_policy == "table_starts"):
        start_codons = start_codons_for(genetic_code, start_codon_policy)

    start_codons = resolve_start_codons(
        start_codons,
        start_codon_policy,
    )

    if stop_codons is None:
        stop_codons = _stop_codons_of(codon_table)

    start_codons = set(
        start_codons
    )

    stop_codons = set(
        stop_codons
    )

    if min_length < 0:

        raise ValueError(
            "min_length must be "
            "greater than or equal to 0."
        )

    reverse_dna = (
        get_reverse_complement(dna)
    )

    strands = [
        ("+", dna),
        ("-", reverse_dna)
    ]

    sequence_length = len(dna)

    orfs = []

    for strand, sequence in strands:

        for frame_offset in range(3):

            frame_number = (
                frame_offset + 1
                if strand == "+"
                else -(frame_offset + 1)
            )

            position = frame_offset

            while (
                position
                <= len(sequence) - 3
            ):

                codon = sequence[
                    position:
                    position + 3
                ]

                start_status = (
                    analyze_codon_identity(
                        codon,
                        start_codons
                    )
                )

                if start_status == "definite":

                    stop_position = (
                        position + 3
                    )

                    stop_codon = None

                    possible_stop_sites = []

                    while (
                        stop_position
                        <= len(sequence) - 3
                    ):

                        current_codon = (
                            sequence[
                                stop_position:
                                stop_position + 3
                            ]
                        )

                        stop_status = (
                            analyze_codon_identity(
                                current_codon,
                                stop_codons
                            )
                        )

                        if (
                            stop_status
                            == "definite"
                        ):

                            stop_codon = (
                                current_codon
                            )

                            break

                        if (
                            stop_status
                            == "possible"
                        ):

                            possible_stop_sites.append({
                                "codon":
                                    current_codon,

                                "position":
                                    stop_position + 1
                            })

                        stop_position += 3

                    # ====================================================
                    # COMPLETE ORF
                    # ====================================================

                    if stop_codon is not None:

                        orf_sequence = sequence[
                            position:
                            stop_position + 3
                        ]

                        if strand == "+":

                            genomic_start = (
                                position
                            )

                            genomic_stop = (
                                stop_position + 2
                            )

                        else:

                            genomic_start = (
                                sequence_length
                                - (
                                    stop_position
                                    + 3
                                )
                            )

                            genomic_stop = (
                                sequence_length
                                - position
                                - 1
                            )

                        if (
                            len(orf_sequence)
                            >= min_length
                        ):

                            ambiguity = (
                                analyze_orf_ambiguity(
                                    orf_sequence, codon_table=codon_table
                                )
                            )

                            codon_count = (
                                len(orf_sequence)
                                // 3
                            )

                            protein_length = (
                                codon_count - 1
                            )

                            start_type = (
                                "Canonical"
                                if codon == "ATG"
                                else "Alternative"
                            )

                            start_context = (
                                analyze_start_context(
                                    sequence,
                                    position
                                )
                            )

                            codon_usage = (
                                analyze_codon_usage(
                                    orf_sequence, codon_table=codon_table
                                )
                            )

                            orf = {
                                "genetic_code":
                                    genetic_code,
                                "strand":
                                    strand,

                                "frame":
                                    frame_number,

                                "start":
                                    genomic_start,

                                "stop":
                                    genomic_stop,

                                "length":
                                    len(orf_sequence),

                                "codon_count":
                                    codon_count,

                                "protein_length":
                                    protein_length,

                                "start_codon":
                                    codon,

                                "start_type":
                                    start_type,

                                "start_status":
                                    start_status,

                                "stop_codon":
                                    stop_codon,

                                "status":
                                    "Complete",

                                "sequence":
                                    orf_sequence,

                                "possible_stop_sites":
                                    possible_stop_sites,

                                "ambiguous_codon_count":
                                    ambiguity[
                                        "ambiguous_codon_count"
                                    ],

                                "ambiguous_codons":
                                    ambiguity[
                                        "ambiguous_codons"
                                    ],

                                "deterministic_ambiguous_count":
                                    ambiguity[
                                        "deterministic_ambiguous_count"
                                    ],

                                "uncertain_codon_count":
                                    ambiguity[
                                        "uncertain_codon_count"
                                    ],

                                "possible_stop_codon_count":
                                    ambiguity[
                                        "possible_stop_codon_count"
                                    ],

                                "translation_safe":
                                    ambiguity[
                                        "translation_safe"
                                    ],

                                "start_context":
                                    start_context,

                                "codon_usage":
                                    codon_usage,

                                "codon_diversity":
                                    calculate_codon_diversity(
                                        orf_sequence, codon_table=codon_table
                                    )
                            }

                            orf["cds_quality"] = (
                                calculate_cds_quality(
                                    orf
                                )
                            )

                            orf["heuristic_score"] = (
                                calculate_orf_heuristic_score(
                                    orf
                                )
                            )

                            orfs.append(orf)

                    # ====================================================
                    # PARTIAL ORF
                    # ====================================================

                    elif include_partial:

                        remaining_length = (
                            len(sequence)
                            - position
                        )

                        usable_length = (
                            remaining_length
                            - (
                                remaining_length
                                % 3
                            )
                        )

                        if usable_length >= 3:

                            orf_sequence = sequence[
                                position:
                                position + usable_length
                            ]

                            if (
                                len(orf_sequence)
                                >= min_length
                            ):

                                last_position = (
                                    position
                                    + usable_length
                                    - 1
                                )

                                if strand == "+":

                                    genomic_start = (
                                        position
                                    )

                                    genomic_stop = (
                                        last_position
                                    )

                                else:

                                    genomic_start = (
                                        sequence_length
                                        - (
                                            position
                                            + usable_length
                                        )
                                    )

                                    genomic_stop = (
                                        sequence_length
                                        - position
                                        - 1
                                    )

                                ambiguity = (
                                    analyze_orf_ambiguity(
                                        orf_sequence, codon_table=codon_table
                                    )
                                )

                                codon_count = (
                                    len(orf_sequence)
                                    // 3
                                )

                                protein_length = (
                                    codon_count
                                )

                                start_type = (
                                    "Canonical"
                                    if codon == "ATG"
                                    else "Alternative"
                                )

                                start_context = (
                                    analyze_start_context(
                                        sequence,
                                        position
                                    )
                                )

                                codon_usage = (
                                    analyze_codon_usage(
                                        orf_sequence, codon_table=codon_table
                                    )
                                )

                                orf = {
                                    "genetic_code":
                                        genetic_code,
                                    "strand":
                                        strand,

                                    "frame":
                                        frame_number,

                                    "start":
                                        genomic_start,

                                    "stop":
                                        genomic_stop,

                                    "length":
                                        len(orf_sequence),

                                    "codon_count":
                                        codon_count,

                                    "protein_length":
                                        protein_length,

                                    "start_codon":
                                        codon,

                                    "start_type":
                                        start_type,

                                    "start_status":
                                        start_status,

                                    "stop_codon":
                                        None,

                                    "status":
                                        "Partial - Missing Stop",

                                    "sequence":
                                        orf_sequence,

                                    "possible_stop_sites":
                                        possible_stop_sites,

                                    "ambiguous_codon_count":
                                        ambiguity[
                                            "ambiguous_codon_count"
                                        ],

                                    "ambiguous_codons":
                                        ambiguity[
                                            "ambiguous_codons"
                                        ],

                                    "deterministic_ambiguous_count":
                                        ambiguity[
                                            "deterministic_ambiguous_count"
                                        ],

                                    "uncertain_codon_count":
                                        ambiguity[
                                            "uncertain_codon_count"
                                        ],

                                    "possible_stop_codon_count":
                                        ambiguity[
                                            "possible_stop_codon_count"
                                        ],

                                    "translation_safe":
                                        ambiguity[
                                            "translation_safe"
                                        ],

                                    "start_context":
                                        start_context,

                                    "codon_usage":
                                        codon_usage,

                                    "codon_diversity":
                                        calculate_codon_diversity(
                                            orf_sequence, codon_table=codon_table
                                        )
                                }

                                orf["cds_quality"] = (
                                    calculate_cds_quality(
                                        orf
                                    )
                                )

                                orf["heuristic_score"] = (
                                    calculate_orf_heuristic_score(
                                        orf
                                    )
                                )

                                orfs.append(orf)

                position += 3

    if str(topology).lower() == "circular":
        orfs = [orf for orf in orfs if orf["status"] == "Complete"]
        orfs.extend(_origin_spanning_orfs(
            dna, min_length=min_length, start_codons=start_codons,
            stop_codons=stop_codons, ranking=ranking, genetic_code=genetic_code))
    elif include_5prime_partial:
        orfs.extend(_five_prime_partial_orfs(
            dna, min_length=min_length, start_codons=start_codons,
            stop_codons=stop_codons, codon_table=codon_table,
            genetic_code=genetic_code, include_partial=include_partial))

    if not nested:
        orfs = _keep_longest_per_stop(orfs)

    # ========================================================
    # RANK ORFS
    # ========================================================

    orfs.sort(
        key=lambda orf: _orf_sort_key(orf, ranking),
        reverse=True
    )

    for rank, orf in enumerate(
        orfs,
        start=1
    ):

        orf["rank"] = rank

    return orfs
