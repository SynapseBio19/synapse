"""Net charge at a pH and the isoelectric point."""

from __future__ import annotations

import math

from analyses._protein.constants import CTERM_PKA, NEGATIVE_PKA, NTERM_PKA, POSITIVE_PKA
from analyses._protein.sequence import get_resolved_protein


def _get_pka_tables(protein):

    positive = POSITIVE_PKA.copy()
    negative = NEGATIVE_PKA.copy()

    if protein:

        first = protein[0]
        last = protein[-1]

        if first in NTERM_PKA:

            positive["Nterm"] = (
                NTERM_PKA[first]
            )

        if last in CTERM_PKA:

            negative["Cterm"] = (
                CTERM_PKA[last]
            )

    return positive, negative


def calculate_charge_at_ph(
    protein,
    pH
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return 0.0

    try:
        pH = float(pH)
    except (TypeError, ValueError):
        raise ValueError("pH must be a numeric value.")

    if not math.isfinite(pH):
        raise ValueError("pH must be finite.")

    if not 0.0 <= pH <= 14.0:
        raise ValueError("pH must be between 0 and 14.")

    positive_pka, negative_pka = (
        _get_pka_tables(
            protein
        )
    )

    positive_charge = 0.0
    negative_charge = 0.0

    # Henderson-Hasselbalch form for basic groups:
    # fraction protonated = 1 / (1 + 10^(pH-pKa))
    for group, pka in positive_pka.items():

        count = (
            1
            if group == "Nterm"
            else protein.count(group)
        )

        if count:
            positive_charge += (
                count
                / (
                    1.0
                    + 10.0 ** (pH - pka)
                )
            )

    # Henderson-Hasselbalch form for acidic groups:
    # fraction deprotonated = 1 / (1 + 10^(pKa-pH))
    for group, pka in negative_pka.items():

        count = (
            1
            if group == "Cterm"
            else protein.count(group)
        )

        if count:
            negative_charge += (
                count
                / (
                    1.0
                    + 10.0 ** (pka - pH)
                )
            )

    return (
        positive_charge
        - negative_charge
    )


def calculate_isoelectric_point(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return None

    # The Bjellqvist pKa model is monotonic in net charge over
    # the pH range 0-14 for the sequence-level Henderson-
    # Hasselbalch model used here. Bisection therefore gives
    # the theoretical pI where calculated net charge is zero.
    low = 0.0
    high = 14.0

    low_charge = calculate_charge_at_ph(
        protein,
        low
    )
    high_charge = calculate_charge_at_ph(
        protein,
        high
    )

    if low_charge < 0 or high_charge > 0:
        return None

    for _ in range(200):

        middle = (
            low + high
        ) / 2.0

        charge = calculate_charge_at_ph(
            protein,
            middle
        )

        if charge > 0:
            low = middle
        else:
            high = middle

    return (
        low + high
    ) / 2.0
