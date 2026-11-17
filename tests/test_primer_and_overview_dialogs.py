"""Single-primer analysis dialog and the overview of a long sequence (Phase 7, D).

The single-primer branch of the primer menu and the profile sections of the
sequence overview (drawn only for sequences of at least 100 bp) were not
run by any test. The expected Tm values are the module's own results for
the same conditions; the comparison with Primer3 is in
`tests/test_primer_thermodynamics.py`.
"""

import builtins
import contextlib
import io
import unittest
from pathlib import Path
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from tests.test_main_session import _Script, run_session
from analyses import primer as P
from ui import primer_menu, style

PRIMER = "AGCGGATAACAATTTCACACAGGA"   # M13 reverse-style 24-mer
FIXTURES = Path(__file__).parent / "fixtures"


def drive(answers, call):
    script = _Script(answers)
    output = io.StringIO()
    with mock.patch.object(builtins, "input", script), contextlib.redirect_stdout(output):
        result = call()
    return output.getvalue(), script.answers, result


class SinglePrimerDialogTests(unittest.TestCase):
    def test_empty_and_invalid_primers_are_errors(self):
        answers = ["3", "",                       # empty
                   "3", "ACGT-12",                # not DNA
                   "3", PRIMER, "9", "0",         # bad model choice, then cancel
                   "3", PRIMER, "1", "abc", "", "", "",   # salt is not a number (all four are read first)
                   "0"]
        output, unused, result = drive(answers, lambda: primer_menu.run_primer_analysis_ui("ACGT" * 60, "demo"))
        self.assertEqual(unused, [])
        self.assertIn("ERROR: Primer sequence is empty.", output)
        self.assertGreaterEqual(output.count("ERROR:"), 3)
        self.assertIn(style.INVALID_SELECTION, output)
        self.assertIsNone(result)

    def test_nearest_neighbor_analysis_with_default_conditions(self):
        answers = ["3", PRIMER, "", "", "", "", ""]   # Enter for the model and all four conditions
        output, unused, result = drive(answers, lambda: primer_menu.run_primer_analysis_ui("ACGT" * 60, "demo"))
        self.assertEqual(unused, [])
        expected = P.analyze_primer_professional(P.normalize_dna_sequence(PRIMER),
                                                 conditions=P.ThermodynamicConfig())
        self.assertIn("PRIMER ANALYSIS", output)
        self.assertIn(f"Nearest-neighbor Tm:   {expected['nearest_neighbor_tm']:.2f} °C", output)
        self.assertIn("Hairpin stem:", output)
        assert isinstance(result, dict)          # the single-primer branch returns its analysis
        self.assertEqual(result["sequence"], PRIMER)

    def test_wallace_analysis_with_given_conditions(self):
        answers = ["3", PRIMER, "2", "100", "2.0", "0.4", "500"]
        output, unused, result = drive(answers, lambda: primer_menu.run_primer_analysis_ui("ACGT" * 60, "demo"))
        self.assertEqual(unused, [])
        self.assertIn("Tm model:", output)
        assert isinstance(result, dict)
        self.assertEqual(result["tm"], P.analyze_primer(PRIMER, tm_method="wallace")["tm"])


class LongSequenceOverviewTests(unittest.TestCase):
    def test_profiles_cpg_and_homopolymers_of_puc19(self):
        answers = [str(FIXTURES / "pUC19_L09137.2.fasta"), "1", "1", "0"]
        output, unused = run_session(answers)
        self.assertEqual(unused, [])
        self.assertIn("GC CONTENT PROFILE", output)
        self.assertIn("SEQUENCE COMPLEXITY PROFILE", output)
        self.assertNotIn("is not available for this sequence", output)
        self.assertIn("Sequence Length:       2686 bp", output)
        self.assertIn("CPG ISLANDS", output)


if __name__ == "__main__":
    unittest.main()
