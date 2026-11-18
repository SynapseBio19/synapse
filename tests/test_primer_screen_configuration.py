"""The primer screening score depends only on its arguments (Phase 7.1, F3).

Until Phase 7.1 ``_v5_candidate_quality`` read the optimum Tm and GC from a
module global (``_V5_ACTIVE_CONFIG``) that the candidate generator set on
every call, so a score depended on the last design that ran. The
configuration is now an argument.
"""

from __future__ import annotations

import dataclasses
import unittest

from analyses import primer as P


def _scored(primer: str, config: P.ProfessionalPrimerConfig) -> P.AdvancedPrimerCandidate:
    relaxed = P._v6_relaxed_generation_config(config.primer)
    candidate = P._make_candidate(primer, "FORWARD", 0, len(primer), relaxed, config.tm_method)
    advanced = P.analyze_primer_candidate_thermodynamics(candidate, conditions=config.thermodynamics)
    return P._v7_apply_selected_tm_model(advanced, candidate, config)


class ScreenConfigurationTests(unittest.TestCase):
    PRIMER = "GTAAAACGACGGCCAGTCA"

    def test_the_optimum_of_the_given_configuration_is_used(self):
        base = P.ProfessionalPrimerConfig()
        shifted = dataclasses.replace(base, primer=dataclasses.replace(
            base.primer, optimal_tm=base.primer.optimal_tm + 3.0, optimal_gc=base.primer.optimal_gc - 5.0))
        item = _scored(self.PRIMER, base)
        tm = float(item.thermodynamics.tm_C)
        gc = float(item.candidate.gc_percent)
        expected = (2.0 * (abs(tm - shifted.primer.optimal_tm) - abs(tm - base.primer.optimal_tm))
                    + (abs(gc - shifted.primer.optimal_gc) - abs(gc - base.primer.optimal_gc)))
        difference = P._v5_candidate_quality(item, base) - P._v5_candidate_quality(item, shifted)
        self.assertAlmostEqual(difference, expected, places=9)

    def test_an_earlier_design_does_not_change_a_later_score(self):
        base = P.ProfessionalPrimerConfig()
        item = _scored(self.PRIMER, base)
        before = P._v5_candidate_quality(item, base)
        other = dataclasses.replace(base, primer=dataclasses.replace(base.primer, optimal_tm=70.0))
        P.design_primers_professional("ACGT" * 60 + "GGCCAATT" * 20, 40, 260, config=other)
        self.assertEqual(P._v5_candidate_quality(item, base), before)
        self.assertEqual(P._v7_screening_score(item.candidate, base), before)

    def test_the_module_keeps_no_active_configuration(self):
        self.assertFalse(hasattr(P, "_V5_ACTIVE_CONFIG"))


if __name__ == "__main__":
    unittest.main()
