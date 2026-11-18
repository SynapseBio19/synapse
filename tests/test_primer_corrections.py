"""Primer-design corrections found on 2026-11-12 (D20, D21).

D20: the candidate screening score (_screening_quality) is documented to
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


# The configuration of the screen; passed explicitly since Phase 7.1 (F3).
CONFIG = P.ProfessionalPrimerConfig()


def _advanced(primer: str) -> P.AdvancedPrimerCandidate:
    config = CONFIG
    relaxed = P._measurement_config(config.primer)
    candidate = P._make_candidate(primer, "FORWARD", 0, len(primer), relaxed, config.tm_method)
    advanced = P.analyze_primer_candidate_thermodynamics(candidate, conditions=config.thermodynamics)
    return P._apply_selected_tm_model(advanced, candidate, config)


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
        # _structure_penalty(8, 3, 8) is the full 30 points.
        self.assertAlmostEqual(
            P._screening_quality(without_stem, CONFIG) - P._screening_quality(advanced, CONFIG), 30.0)

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
        expected = P._structure_penalty(stem, 4, 9)
        self.assertGreater(expected, 0.0)
        self.assertAlmostEqual(
            P._screening_quality(without_stem, CONFIG) - P._screening_quality(advanced, CONFIG), expected)


class PairAssessmentExportTests(unittest.TestCase):
    """D23: the ranking's risk flags and warnings were attributes outside the
    dataclass fields, so the JSON export (which serializes fields) lost them."""

    def test_the_pair_assessment_reaches_the_export(self):
        import json
        from core import api
        from core.report import canonical_json
        from tests._support import deterministic_dna
        result = api.design_primers(deterministic_dna(500))
        exported = json.loads(canonical_json(result.data))["ranked_pairs"]
        self.assertTrue(exported)
        for advanced, record in zip(result.data.ranked_pairs, exported):
            for name in ("screening_score", "high_risk", "critical_risk",
                         "structural_warnings", "specificity_screened"):
                self.assertEqual(record[name], getattr(advanced, name), name)
            self.assertAlmostEqual(record["model_tm_delta"], advanced.model_tm_delta)


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
