"""P2-6: the cloning planner computes each enzyme analysis once per search.

The planner evaluates every enzyme combination over the same two molecules.
Before the search cache each of the ~63,500 evaluations for 22 enzymes
re-validated both sequences and re-scanned them for every enzyme (about
258,000 scans, 21.7 s for a 912 bp vector and a 252 bp insert). The cache
lives only for one search, returns independent copies, and must not change
a single result.
"""

from __future__ import annotations

import unittest
from unittest import mock

from tests._support import deterministic_dna

VECTOR = deterministic_dna(400) + "GAATTC" + deterministic_dna(300)[::-1] + "GGATCC" + deterministic_dna(200)
INSERT = "GAATTC" + deterministic_dna(240) + "GGATCC"
ENZYMES = ["EcoRI", "BamHI", "HindIII", "XhoI", "PstI"]


class PlannerCacheTests(unittest.TestCase):
    def test_each_analysis_is_computed_once_per_search(self):
        from analyses import cloning as C
        calls = []
        original = C._analyze_enzyme

        def counting(sequence, enzyme, topology):
            calls.append((len(sequence), enzyme, topology))
            return original(sequence, enzyme, topology)

        with mock.patch.object(C, "_analyze_enzyme", counting):
            result = C.automatic_cloning_planner(VECTOR, INSERT, enzyme_names=ENZYMES, top_n=3)
        self.assertGreaterEqual(result["evaluated_count"], 30)    # each needs four analyses
        self.assertEqual(len(calls), len(set(calls)))            # no analysis twice
        self.assertLessEqual(len(calls), 2 * len(ENZYMES))       # two molecules x enzymes

    def test_results_equal_an_uncached_search(self):
        from analyses import cloning as C
        for topology in ("linear", "circular"):
            with self.subTest(vector_topology=topology):
                cached = C.automatic_cloning_planner(VECTOR, INSERT, enzyme_names=ENZYMES, top_n=50,
                                                     vector_topology=topology)
                self.assertIsNone(C._STRATEGY_SEARCH_CACHE)     # the cache ended with the search
                with mock.patch.object(C, "_cached", lambda key, compute: compute()):
                    uncached = C.automatic_cloning_planner(VECTOR, INSERT, enzyme_names=ENZYMES, top_n=50,
                                                           vector_topology=topology)
                self.assertEqual(cached, uncached)

    def test_results_share_no_objects(self):
        from analyses import cloning as C
        with C._strategy_search_cache():
            first = C.analyze_enzyme(VECTOR, "EcoRI")
            first["sites"].append(-1)
            first["site_matches"].clear()
            second = C.analyze_enzyme(VECTOR, "EcoRI")
        self.assertNotIn(-1, second["sites"])
        self.assertTrue(second["site_matches"])
        self.assertEqual(second, C.analyze_enzyme(VECTOR, "EcoRI"))   # same as uncached


if __name__ == "__main__":
    unittest.main()
