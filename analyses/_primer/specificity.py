"""Specificity screen against caller-supplied reference sequences and the
binding sites of a pair. No genome-wide claim is made."""

from __future__ import annotations

from typing import Optional, Sequence

from analyses._primer.sequence import normalize_dna_sequence, reverse_complement
from analyses._primer.models import SpecificityHit


# ============================================================================
# SPECIFICITY SCREENING
# ============================================================================

def screen_primer_specificity(
    primer: str,
    reference_sequences: Optional[dict[str, str] | Sequence[tuple[str, str]]] = None,
    max_mismatches: int = 2,
    three_prime_window: int = 5,
) -> list[SpecificityHit]:
    """Screen a primer against supplied references; this is not genome-wide BLAST."""
    primer = normalize_dna_sequence(primer)
    if max_mismatches < 0:
        raise ValueError("max_mismatches must be non-negative.")
    if three_prime_window <= 0:
        raise ValueError("three_prime_window must be positive.")
    if not reference_sequences:
        return []
    references = list(reference_sequences.items()) if isinstance(reference_sequences, dict) else list(reference_sequences)
    hits = []
    for name, raw_reference in references:
        reference = normalize_dna_sequence(raw_reference)
        orientations = (
            ("FORWARD_MATCH", primer),
            ("REVERSE_COMPLEMENT_MATCH", reverse_complement(primer)),
        )
        for orientation, query in orientations:
            if len(query) > len(reference):
                continue
            for start in range(len(reference) - len(query) + 1):
                window = reference[start:start + len(query)]
                mismatches = sum(a != b for a, b in zip(query, window))
                if mismatches > max_mismatches:
                    continue
                terminal_length = min(three_prime_window, len(query))
                three_prime_mismatches = sum(
                    a != b for a, b in zip(query[-terminal_length:], window[-terminal_length:])
                )
                hits.append(SpecificityHit(
                    sequence_name=name,
                    start=start,
                    end=start + len(query),
                    mismatches=mismatches,
                    three_prime_mismatches=three_prime_mismatches,
                    exact=(mismatches == 0),
                    orientation=orientation,
                ))
    hits.sort(key=lambda hit: (hit.mismatches, hit.three_prime_mismatches))
    return hits


# ============================================================================
# SYNAPSE — PRIMER DESIGN UI
# ============================================================================
#
# Module-owned UI.
# The primer engine above is intentionally preserved.
# ============================================================================


def primer_pair_binding_sites(pair, three_prime_window: int = 5) -> list[dict]:
    """The two binding sites of a primer pair, 1-based inclusive (audit P2-4).

    Accepts a PrimerPair or an AdvancedPrimerPair. A primer's start/end are
    0-based half-open on the template; the forward primer's 3' end is the
    right end of its site, the reverse primer's 3' end the left end (it binds
    the top strand reverse-complemented). three_prime_window bases at the 3'
    end form the critical region, where a mismatch most affects extension.
    """
    pair = getattr(pair, "pair", pair)
    window = max(int(three_prime_window), 1)
    forward, reverse = pair.forward, pair.reverse
    f_start, f_end = forward.start + 1, forward.end
    r_start, r_end = reverse.start + 1, reverse.end
    return [
        {"name": "forward primer", "kind": "primer", "start": f_start, "end": f_end,
         "critical": (max(f_start, f_end - window + 1), f_end), "critical_label": "3' end"},
        {"name": "reverse primer", "kind": "primer", "start": r_start, "end": r_end,
         "critical": (r_start, min(r_end, r_start + window - 1)), "critical_label": "3' end"},
    ]
