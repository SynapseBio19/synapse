"""Stage 5: unresolved residues ("?") must not shift position-dependent results."""

from __future__ import annotations

import unittest

from tests import _support  # noqa: F401

UBIQUITIN = ("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
AMBIGUOUS = UBIQUITIN[:30] + "?" + UBIQUITIN[31:]


class PositionalProteinTests(unittest.TestCase):
    def test_positional_helper_keeps_question_marks(self):
        from analyses import protein as P
        # Every position is kept: the internal "*" and the ambiguity code X
        # become "?" (Phase 6, D14-F6; previously both were deleted).
        self.assertEqual(P.get_positional_protein("MK?L*x"), "MK?L??")
        self.assertEqual(P.get_positional_protein("MK?L*"), "MK?L")  # terminal stop
        self.assertEqual(P.get_resolved_protein("MK?L"), "MKL")

    def test_window_positions_are_preserved(self):
        from analyses import protein as P
        clean = {w["start"]: round(w["score"], 9) for w in P.calculate_hydropathy_profile(UBIQUITIN, 9)}
        ambiguous = {w["start"]: round(w["score"], 9) for w in P.calculate_hydropathy_profile(AMBIGUOUS, 9)}
        self.assertEqual(len(clean), 68)
        self.assertEqual(len(ambiguous), 68 - 9)  # only the 9 windows covering the "?" are dropped
        for start, score in ambiguous.items():
            self.assertEqual(score, clean[start], f"window at {start} shifted")

    def test_flexibility_and_signal_peptide_skip_unresolved_windows(self):
        from analyses import protein as P
        flexible = P.calculate_flexibility_profile(AMBIGUOUS, 9)
        self.assertTrue(all("?" not in w["sequence"] for w in flexible))
        self.assertTrue(P.predict_signal_peptide("MKLLLLAAALLLLACSSA" + "K" * 20)["predicted"])
        self.assertFalse(P.predict_signal_peptide("M?LLLLAAAGGCSS" + "K" * 20)["predicted"])

    def test_instability_index_ignores_unresolved_dipeptides(self):
        from analyses import protein as P
        index = P.calculate_instability_index(UBIQUITIN)
        assert index is not None
        self.assertAlmostEqual(index, 36.055263, places=5)
        self.assertIsNotNone(P.calculate_instability_index(AMBIGUOUS))
        self.assertIsNone(P.calculate_instability_index("??"))

    def test_translation_with_ambiguous_codon_reports_real_positions(self):
        from analyses import protein as P
        analysis = P.analyze_complete_protein("", rna="AUG" + "RAA" + "CUG" * 20 + "UAA",
                                              initiator_start=True)
        self.assertEqual(analysis["protein"][:3], "M?L")
        self.assertEqual(analysis["transmembrane_domains"][0]["start"], 3)
        self.assertTrue(all("?" not in w["sequence"] for w in analysis["hydropathy"]))


if __name__ == "__main__":
    unittest.main()
