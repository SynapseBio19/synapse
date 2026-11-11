"""P2-5: G-U wobble pairs in the RNA hairpin screen.

G-U is a regular pair in RNA helices, so a stem such as GGUGC/GUAUC (three
Watson-Crick and two wobble pairs) is a real stem-loop candidate that a
Watson-Crick-only screen misses. Counting G-U like G-C without free
energies, however, multiplies the candidates (58 -> 260 on the golden-master
reference), so the screen keeps Watson-Crick pairing as its rule and reports
the G-U count beside it.
"""

from __future__ import annotations

import unittest

from tests import _support  # noqa: F401

WOBBLE_STEM = "GGUGCAAAAGUAUC"     # GGUGC with CUAUG read backwards: G-C, G-U, U-A, G-U, C-G


class WobblePairTests(unittest.TestCase):
    def test_a_wobble_stem_is_found_only_when_g_u_counts(self):
        from analyses.rna import find_rna_hairpins
        self.assertEqual(find_rna_hairpins(WOBBLE_STEM), [])
        (hairpin,) = [h for h in find_rna_hairpins(WOBBLE_STEM, allow_gu=True)
                      if (h["start"], h["end"], h["stem_length"]) == (1, 14, 5)]
        self.assertEqual(hairpin["loop_length"], 4)
        self.assertEqual((hairpin["pairing_fraction"], hairpin["wobble_pairs"]), (1.0, 2))

    def test_watson_crick_stems_carry_no_wobble_pairs(self):
        from analyses.rna import find_rna_hairpins
        for hairpin in find_rna_hairpins("GGGCAAAAGCCC"):
            self.assertEqual(hairpin["wobble_pairs"], 0)

    def test_summary_reports_both_counts_and_keeps_the_screen(self):
        from analyses.rna import summarize_rna_structure
        from tests.golden_scenarios import REFERENCE
        summary = summarize_rna_structure(REFERENCE.replace("T", "U"))
        self.assertEqual(summary["hairpin_count"], 58)            # unchanged screen
        self.assertEqual(summary["hairpin_count_with_gu"], 260)   # with wobble pairs
        self.assertEqual(summarize_rna_structure(WOBBLE_STEM)["hairpin_count"], 0)
        self.assertGreater(summarize_rna_structure(WOBBLE_STEM)["hairpin_count_with_gu"], 0)


if __name__ == "__main__":
    unittest.main()
