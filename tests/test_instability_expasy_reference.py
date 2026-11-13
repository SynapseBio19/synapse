"""Instability-index weights compared with ExPASy ProtParam (measured 2026-11-13).

ProtParam computes the instability index as II = (10 / L) * sum of the
weights of the L-1 dipeptides (Guruprasad et al. 1990). Querying it with
the 20 probes ``X*5`` and the 380 probes ``X*4 + Y`` therefore determines
every one of the 400 dipeptide weights it uses:

    II(X*5)     = (10/5) * 4 * w[XX]                 -> w[XX]
    II(X*4 + Y) = (10/5) * (3 * w[XX] + w[XY])       -> w[XY]

The measured II values are stored in
``tests/fixtures/expasy/protparam_instability_probes.json``; the weights
are recovered here, so nothing of ProtParam's own table is copied.

Result: 399 of 400 weights agree with Synapse's table. The single
difference is KQ, where ProtParam uses 24.68 and Synapse 24.64 (the value
Biopython's ``ProtParamData.DIWV`` also carries). Which of the two is the
value of the 1990 paper cannot be decided without the paper itself, so the
table is left unchanged and the difference is recorded (docs/METHODS.md).
Its effect is bounded: one KQ dipeptide shifts the index of a 100-residue
protein by 0.004.
"""

import hashlib
import json
import unittest
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from analyses.protein import INSTABILITY_WEIGHTS, calculate_instability_index

FIXTURE = Path(__file__).parent / "fixtures" / "expasy" / "protparam_instability_probes.json"
FIXTURE_SHA256 = "146cc390200676e3a666d5bb6c9601117ac49a2e2b6dd7a83a70367aa79e6965"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
KNOWN_DIFFERENCES = {"KQ": (24.68, 24.64)}  # ProtParam, Synapse/Biopython
TOLERANCE = 0.02  # ProtParam prints the index with two decimals.


def measured():
    return json.loads(FIXTURE.read_text())["instability_index"]


def recovered_weights(values):
    """The 400 weights ProtParam uses, from the probe measurements."""
    weights = {}
    for first in AMINO_ACIDS:
        weights[first + first] = values[first * 5] * 5 / 40
    for first in AMINO_ACIDS:
        for second in AMINO_ACIDS:
            if first == second:
                continue
            index = values[first * 4 + second]
            weights[first + second] = index * 5 / 10 - 3 * weights[first + first]
    return weights


class ExpasyInstabilityWeightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = measured()
        cls.weights = recovered_weights(cls.values)

    def test_fixture_is_unchanged(self):
        self.assertEqual(hashlib.sha256(FIXTURE.read_bytes()).hexdigest(), FIXTURE_SHA256)
        self.assertEqual(len(self.values), 400)

    def test_every_weight_agrees_except_the_recorded_difference(self):
        differences = {}
        for dipeptide, weight in self.weights.items():
            synapse = INSTABILITY_WEIGHTS[dipeptide[0]][dipeptide[1]]
            if abs(weight - synapse) > TOLERANCE:
                differences[dipeptide] = (round(weight, 2), synapse)
        self.assertEqual(differences, KNOWN_DIFFERENCES)

    def test_probe_indices_are_reproduced(self):
        """Synapse's own index for each probe equals ProtParam's, except where KQ occurs."""
        for probe, expasy in self.values.items():
            with self.subTest(probe=probe):
                difference = abs(calculate_instability_index(probe) - expasy)
                # One KQ in a 5-residue probe shifts the index by 10/5 * 0.04 = 0.08.
                self.assertLessEqual(difference, 0.09 if "KQ" in probe else 0.005)


if __name__ == "__main__":
    unittest.main()
