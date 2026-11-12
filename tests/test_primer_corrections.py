"""Primer-design corrections found on 2026-11-12 (D20, D21).

D20: the candidate screening score (_v5_candidate_quality) is documented to
penalize hairpin and self-dimer stems, but it read the stems from the
thermodynamic profile, which has no such fields. The penalty was always
zero. Both tests below fail on the code before the correction.

D21: ProfessionalPrimerConfig accepts tm_method "nearest-neighbor", but
candidate generation raised "Unknown Tm method" for it.
"""

from __future__ import annotations

import dataclasses
import unittest

from analyses import primer as P


def _advanced(primer: str) -> P.AdvancedPrimerCandidate:
    config = P.ProfessionalPrimerConfig()
    P._V5_ACTIVE_CONFIG = config
    relaxed = P._v6_relaxed_generation_config(config.primer)
    candidate = P._make_candidate(primer, "FORWARD", 0, len(primer), relaxed, config.tm_method)
    advanced = P.analyze_primer_candidate_thermodynamics(candidate, conditions=config.thermodynamics)
    return P._v7_apply_selected_tm_model(advanced, candidate, config)


class StructuralPenaltyTests(unittest.TestCase):
    def test_a_hairpin_stem_lowers_the_screening_score(self):
        # GCGCATGA ... TCATGCGC closes an 8-bp stem around a 3-nt loop.
        advanced = _advanced("GCGCATGAAATTCATGCGC")
        self.assertEqual(advanced.candidate.hairpin.max_stem, 8)
        without_stem = dataclasses.replace(
            advanced,
            candidate=dataclasses.replace(
                advanced.candidate,
                hairpin=dataclasses.replace(advanced.candidate.hairpin, max_stem=0),
            ),
        )
        # _v5_structure_penalty(8, 3, 8) is the full 30 points.
        self.assertAlmostEqual(
            P._v5_candidate_quality(without_stem) - P._v5_candidate_quality(advanced), 30.0)

    def test_a_self_dimer_stem_lowers_the_screening_score(self):
        advanced = _advanced("ATGAATTCATGCAGGTCAAG")
        stem = advanced.candidate.self_dimer.max_stem
        self.assertGreaterEqual(stem, 5)
        without_stem = dataclasses.replace(
            advanced,
            candidate=dataclasses.replace(
                advanced.candidate,
                self_dimer=dataclasses.replace(advanced.candidate.self_dimer, max_stem=0),
            ),
        )
        expected = P._v5_structure_penalty(stem, 4, 9)
        self.assertGreater(expected, 0.0)
        self.assertAlmostEqual(
            P._v5_candidate_quality(without_stem) - P._v5_candidate_quality(advanced), expected)


class TmMethodNameTests(unittest.TestCase):
    def test_nearest_neighbor_with_a_hyphen_designs_like_the_underscore_name(self):
        template = "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA"
        hyphen = P.design_primers_professional(
            template, config=P.ProfessionalPrimerConfig(tm_method="nearest-neighbor"))
        underscore = P.design_primers_professional(
            template, config=P.ProfessionalPrimerConfig(tm_method="nearest_neighbor"))
        self.assertTrue(hyphen.ranked_pairs)
        self.assertEqual(
            [(p.pair.forward.sequence, p.pair.reverse.sequence, p.pair.score) for p in hyphen.ranked_pairs],
            [(p.pair.forward.sequence, p.pair.reverse.sequence, p.pair.score) for p in underscore.ranked_pairs],
        )

    def test_calculate_tm_accepts_the_hyphenated_name(self):
        self.assertEqual(P.calculate_tm("ACGTACGTACGTACGTAC", "nearest-neighbor"),
                         P.calculate_tm("ACGTACGTACGTACGTAC", "nearest_neighbor"))


if __name__ == "__main__":
    unittest.main()
