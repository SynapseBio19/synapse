"""Results of the BLAST and primer modules stay identical (Phase 7.1).

See ``tests/characterization.py`` for what is recorded and why.
"""

from __future__ import annotations

import json
import math
import os
import unittest

from tests import characterization as C


# Floats are compared exactly, except for a relative difference of at most
# 1e-12. CPython 3.10.11 and 3.14.6 on Windows return results of math.log2
# that differ in the last digit: sequence_complexity of the "deletion"
# variant is 98.82341858018908 on 3.10 and 98.82341858018907 on 3.14, with
# the unchanged code of 727208a on both (checked 2026-11-18). That is the
# interpreter's math library, not a change of Synapse. Since Python 3.12,
# sum() of floats is compensated (more exact), so sums that should be 0
# differ around 1e-17 between 3.10 and 3.14 (protein hydropathy windows,
# found 2026-11-19); ABS_TOL covers values near zero.
REL_TOL = 1e-12
ABS_TOL = 1e-12


def _as_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _first_difference(actual, expected, path):
    """None if equal within REL_TOL for floats, else the path of the first difference."""
    if actual == expected:
        return None
    if isinstance(actual, dict) and isinstance(expected, dict):
        if sorted(actual) != sorted(expected):
            return f"{path}: keys differ"
        for key in expected:
            found = _first_difference(actual[key], expected[key], f"{path}/{key}")
            if found:
                return found
        return None
    if isinstance(actual, list) and isinstance(expected, list):
        if len(actual) != len(expected):
            return f"{path}: lengths differ"
        for index, (a, e) in enumerate(zip(actual, expected)):
            found = _first_difference(a, e, f"{path}[{index}]")
            if found:
                return found
        return None
    if isinstance(actual, (int, float)) and isinstance(expected, (int, float))             and not isinstance(actual, bool) and not isinstance(expected, bool):
        return None if math.isclose(actual, expected, rel_tol=REL_TOL, abs_tol=ABS_TOL) else f"{path}: {actual!r} != {expected!r}"
    if isinstance(actual, str) and isinstance(expected, str):
        a, e = _as_float(actual), _as_float(expected)
        if a is not None and e is not None:
            return _first_difference(a, e, path)
        if actual[:1] in "{[" and expected[:1] in "{[":
            try:
                return _first_difference(json.loads(actual), json.loads(expected), path + "<json>")
            except ValueError:
                pass
    return f"{path}: {str(actual)[:120]!r} != {str(expected)[:120]!r}"


class CharacterizationSnapshotTests(unittest.TestCase):
    maxDiff = 4000

    @classmethod
    def setUpClass(cls):
        update = os.environ.get("SYNAPSE_UPDATE_CHARACTERIZATION", "")
        if update:
            C.write_snapshot("" if update == "1" else update)

    def _compare(self, section: str, actual: dict) -> None:
        expected = json.loads(C.snapshot_for(section).read_text(encoding="utf-8"))[section]
        # Round-trip through JSON so tuples and lists compare alike.
        actual = json.loads(json.dumps(actual, sort_keys=True))
        self.assertEqual(sorted(actual), sorted(expected))
        for key in expected:
            with self.subTest(section=section, scenario=key):
                if actual[key] != expected[key]:
                    difference = _first_difference(actual[key], expected[key], key)
                    self.assertIsNone(difference, difference)

    def test_primer_results_are_unchanged(self):
        self._compare("primer", C.primer_scenarios())

    def test_blast_results_are_unchanged(self):
        self._compare("blast", C.blast_scenarios())

    def test_mutation_results_are_unchanged(self):
        self._compare("mutation", C.mutation_scenarios())

    def test_cloning_results_are_unchanged(self):
        self._compare("cloning", C.cloning_scenarios())

    def test_protein_results_are_unchanged(self):
        self._compare("protein", C.protein_scenarios())

    def test_variant_results_are_unchanged(self):
        self._compare("variant", C.variant_scenarios())

    def test_sequence_results_are_unchanged(self):
        self._compare("sequence", C.sequence_scenarios())

    def test_orf_results_are_unchanged(self):
        self._compare("orf", C.orf_scenarios())

    def test_restriction_results_are_unchanged(self):
        self._compare("restriction", C.restriction_scenarios())

    def test_rna_results_are_unchanged(self):
        self._compare("rna", C.rna_scenarios())


if __name__ == "__main__":
    unittest.main()
