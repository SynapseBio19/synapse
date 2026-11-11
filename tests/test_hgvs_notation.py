"""Predicted protein consequences in HGVS style (den Dunnen et al. 2016).

A consequence inferred from DNA is a prediction, so HGVS writes it in
parentheses. One-letter codes throughout, "*" for a stop codon (one-letter
style; "Ter" belongs to the three-letter style), "=" for no change, "?" where
the effect is unknown. Before 2026-11-11 the helper wrote p.R2C without
parentheses and mixed "Ter" into one-letter notation, no report used it, and
the variant bridge stored the effect class ("Missense") as protein_change.
"""

from __future__ import annotations

import unittest

from tests import _support  # noqa: F401

REFERENCE = "ATGCGTTGGAAATAA"      # M R W K *


class ProteinNotationTests(unittest.TestCase):
    CASES = {
        "ATGTGTTGGAAATAA": ("Missense", "p.(R2C)"),
        "ATGCGCTGGAAATAA": ("Synonymous", "p.(R2=)"),
        "ATGCGTTGAAAATAA": ("Nonsense", "p.(W3*)"),
        "GTGCGTTGGAAATAA": ("Start-loss", "p.(M1?)"),
        "ATGCGTTGGAAACAA": ("Stop-loss", "p.(*5Qext*?)"),
    }

    def test_each_consequence_class(self):
        from analyses import mutation as M
        for variant, (effect, notation) in self.CASES.items():
            with self.subTest(variant=variant):
                (mutation,) = M.compare_sequences(REFERENCE, variant).mutations
                self.assertEqual((mutation.effect, M.mutation_hgvs_p(mutation)), (effect, notation))

    def test_no_notation_without_a_codon(self):
        from analyses import mutation as M
        (insertion,) = M.compare_sequences(REFERENCE, "ATGCGTAAATGGAAATAA").mutations
        self.assertIsNone(M.mutation_hgvs_p(insertion))

    def test_report_line_names_the_predicted_change(self):
        from analyses import mutation as M
        (mutation,) = M.compare_sequences(REFERENCE, "ATGTGTTGGAAATAA").mutations
        self.assertTrue(M.format_mutation(mutation).endswith("4C>T  [Missense p.(R2C)]"))

    def test_variant_bridge_stores_the_notation_not_the_class(self):
        from analyses import variant_analysis as V
        variant = V.annotate_from_mutation_analysis(
            V.Variant("c", 4, "A", "G"), reference_sequence="ATGAAATAA", alternative_sequence="ATGGAATAA")
        self.assertEqual((variant.consequence, variant.protein_change), ("Missense", "p.(K2E)"))


if __name__ == "__main__":
    unittest.main()
