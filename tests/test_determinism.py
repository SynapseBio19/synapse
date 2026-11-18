"""Results do not depend on Python's hash seed (2026-11-19).

Python randomizes the iteration order of sets per process (PYTHONHASHSEED).
``analyses.sequence.calculate_sequence_complexity`` summed the entropy terms
in the order of the set {"A", "T", "C", "G"}, so the last digit of the
complexity, and with it the reported position of the complexity maximum,
changed from run to run. The same check over all characterization scenarios
(blast, primer, mutation, cloning, protein, variant, sequence) found no
other case. This test runs a sequence analysis under four hash seeds and
requires identical output.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import unittest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROBE = """
import json
from analyses import sequence as S
dna = "CGCGGCGCGTACGCGCGGGCGCCCGCGATCG" * 12 + "AAAAAAAATTTTTGGGGGGCCCC" * 5 + "ACGT" * 40
result = S.analyze_sequence(dna, profile_window=50, profile_step=5)
print(json.dumps(result, sort_keys=True, default=repr))
"""


class HashSeedDeterminismTests(unittest.TestCase):
    def test_sequence_analysis_is_identical_under_any_hash_seed(self):
        outputs = set()
        for seed in ("1", "2", "3", "4"):
            environment = dict(os.environ, PYTHONHASHSEED=seed)
            completed = subprocess.run([sys.executable, "-c", PROBE], cwd=PROJECT_ROOT, env=environment,
                                       capture_output=True, text=True, check=True)
            outputs.add(completed.stdout)
        self.assertEqual(len(outputs), 1, "the result depends on the hash seed")
        self.assertIn("complexity_profile_statistics", json.loads(outputs.pop()))


if __name__ == "__main__":
    unittest.main()
