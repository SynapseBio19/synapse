"""Component tests for the primer modules (Phase 7.1).

After the split of ``analyses/primer.py`` the least-reached parts were the
rejection of invalid reaction conditions (``thermodynamics``), the
specificity screen against supplied references (``specificity``) and the
failure reasons and warnings of the professional workflow (``workflow``).
These tests pin that documented behavior; none of it changes a result.
"""

from __future__ import annotations

import dataclasses
import unittest

from analyses import primer as P

PRIMER = "GTAAAACGACGGCCAGT"


class ConditionValidationTests(unittest.TestCase):
    def test_invalid_conditions_are_rejected_with_their_reason(self):
        base = P.ThermodynamicConfig()
        cases = [
            ({"sodium_mM": 0.0}, "sodium_mM must be greater than zero."),
            ({"magnesium_mM": -0.1}, "magnesium_mM cannot be negative."),
            ({"dntp_mM": -0.1}, "dntp_mM cannot be negative."),
            ({"primer_nM": 0.0}, "primer_nM must be greater than zero."),
            ({"temperature_reference_C": -300.0}, "temperature_reference_C is below absolute zero."),
            ({"dna_concentration_mode": "mixed"}, "dna_concentration_mode must be 'duplex' or 'strand'."),
        ]
        for change, message in cases:
            with self.subTest(change=change):
                with self.assertRaises(ValueError) as raised:
                    dataclasses.replace(base, **change).validate()
                self.assertEqual(str(raised.exception), message)
        base.validate()

    def test_ionic_strength_inputs(self):
        for kwargs, message in (({"sodium_mM": 0.0}, "sodium_mM must be greater than zero."),
                                ({"magnesium_mM": -1.0}, "magnesium_mM cannot be negative."),
                                ({"dntp_mM": -1.0}, "dntp_mM cannot be negative."),
                                ({"salt_model": "debye"}, "Unknown salt model.")):
            arguments = {"sodium_mM": 50.0, "magnesium_mM": 1.5, "dntp_mM": 0.2, **kwargs}
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError) as raised:
                    P.calculate_effective_ionic_strength(**arguments)
                self.assertEqual(str(raised.exception), message)
        # "na_only" ignores Mg2+ and dNTPs by definition.
        self.assertEqual(P.calculate_effective_ionic_strength(50.0, 3.0, 0.8, salt_model="na_only"), 50.0)
        self.assertGreater(P.calculate_effective_ionic_strength(50.0, 3.0, 0.2), 50.0)


class SpecificityScreenTests(unittest.TestCase):
    def test_arguments_are_checked(self):
        with self.assertRaisesRegex(ValueError, "max_mismatches must be non-negative."):
            P.screen_primer_specificity(PRIMER, {"r": "ACGT"}, max_mismatches=-1)
        with self.assertRaisesRegex(ValueError, "three_prime_window must be positive."):
            P.screen_primer_specificity(PRIMER, {"r": "ACGT"}, three_prime_window=0)
        self.assertEqual(P.screen_primer_specificity(PRIMER, None), [])
        self.assertEqual(P.screen_primer_specificity(PRIMER, {"short": "ACGTACGT"}), [])

    def test_both_orientations_and_both_input_forms(self):
        flank = "TTTTTTTTTT"
        forward = flank + PRIMER + flank
        reverse = flank + P.reverse_complement(PRIMER) + flank
        hits = P.screen_primer_specificity(PRIMER, [("fwd", forward), ("rev", reverse)])
        found = {(hit.sequence_name, hit.orientation, hit.start, hit.mismatches) for hit in hits}
        self.assertIn(("fwd", "FORWARD_MATCH", 10, 0), found)
        self.assertIn(("rev", "REVERSE_COMPLEMENT_MATCH", 10, 0), found)
        self.assertEqual(P.screen_primer_specificity(PRIMER, {"fwd": forward, "rev": reverse}), hits)


class ProfessionalWorkflowMessageTests(unittest.TestCase):
    TEMPLATE = ("TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATG"
                "CCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAG"
                "ATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCC"
                "ATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTG")

    def test_without_references_the_report_says_so(self):
        report = P.design_primers_professional(self.TEMPLATE, 80, 260)
        self.assertIn("Specificity was not genome-wide screened; no reference database was supplied.",
                      report.warnings)

    def test_references_are_reported_with_their_hits(self):
        report = P.design_primers_professional(self.TEMPLATE, 80, 260, reference_sequences={"self": self.TEMPLATE})
        self.assertTrue(report.ranked_pairs)
        self.assertTrue(any(w.startswith("Specificity screening found ") for w in report.warnings), report.warnings)
        unrelated = {"other": "AC" * 200}
        report = P.design_primers_professional(self.TEMPLATE, 80, 260, reference_sequences=unrelated)
        self.assertIn("Specificity screening was performed against the supplied reference sequences; "
                      "no significant hits were detected in the reported pairs.", report.warnings)

    def test_a_template_without_a_pair_explains_why(self):
        report = P.design_primers_professional("A" * 60, 10, 50)
        self.assertEqual(report.ranked_pairs, [])
        self.assertEqual(report.failure_reasons,
                         ["No primer pair satisfies the configured amplicon and compatibility constraints."])

    def test_a_homopolymer_pair_is_reported_as_critical(self):
        # The professional screen measures every window without hard filters
        # (see the V6 section comment in analyses/_primer/screening.py) and
        # flags risk in the ranking; it does not hide the best available pair.
        report = P.design_primers_professional("A" * 400, 100, 300)
        self.assertEqual(len(report.ranked_pairs), 1)
        pair = report.ranked_pairs[0]
        self.assertTrue(pair.high_risk)
        self.assertTrue(pair.critical_risk)
        self.assertTrue(any(w.startswith("heterodimer risk") for w in pair.structural_warnings))

if __name__ == "__main__":
    unittest.main()
