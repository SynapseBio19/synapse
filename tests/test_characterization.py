"""Results of the BLAST and primer modules stay identical (Phase 7.1).

See ``tests/characterization.py`` for what is recorded and why.
"""

from __future__ import annotations

import json
import os
import unittest

from tests import characterization as C


class CharacterizationSnapshotTests(unittest.TestCase):
    maxDiff = 4000

    @classmethod
    def setUpClass(cls):
        if os.environ.get("SYNAPSE_UPDATE_CHARACTERIZATION") == "1":
            C.write_snapshot()
        cls.expected = json.loads(C.SNAPSHOT.read_text(encoding="utf-8"))

    def _compare(self, section: str, actual: dict) -> None:
        expected = self.expected[section]
        # Round-trip through JSON so tuples and lists compare alike.
        actual = json.loads(json.dumps(actual, sort_keys=True))
        self.assertEqual(sorted(actual), sorted(expected))
        for key in expected:
            with self.subTest(section=section, scenario=key):
                self.assertEqual(actual[key], expected[key])

    def test_primer_results_are_unchanged(self):
        self._compare("primer", C.primer_scenarios())

    def test_blast_results_are_unchanged(self):
        self._compare("blast", C.blast_scenarios())


if __name__ == "__main__":
    unittest.main()
