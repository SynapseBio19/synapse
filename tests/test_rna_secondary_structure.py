"""Optional energy-based RNA structure (P2-5) and its comparison with the hairpin screen.

ViennaRNA is an optional reference package (docs/TOOL_AND_DATA_LEDGER.md).
Tests that need it skip with the reason when it cannot be loaded; the
behavior without it is tested by hiding the module.
"""

import io
import random
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import rna_structure
from analyses.rna import find_rna_hairpins, transcribe_rna

STATUS = rna_structure.viennarna_status()
NEEDS_VIENNARNA = unittest.skipUnless(STATUS["available"], STATUS["reason"] or "")
PUC19 = Path(__file__).parent / "fixtures" / "pUC19_L09137.2.fasta"


class HairpinParserTests(unittest.TestCase):
    def test_single_hairpin(self):
        self.assertEqual(rna_structure.mfe_hairpins("..((((....))))..."),
                         [{"start": 3, "end": 14, "stem_length": 4, "loop_start": 7,
                           "loop_end": 10, "loop_length": 4}])

    def test_bulge_ends_the_stem_and_two_hairpins_are_found(self):
        hairpins = rna_structure.mfe_hairpins("((.(((...)))))..((...))")
        self.assertEqual([(h["start"], h["end"], h["stem_length"]) for h in hairpins],
                         [(4, 12, 3), (17, 23, 2)])

    def test_unpaired_structure_has_no_hairpin(self):
        self.assertEqual(rna_structure.mfe_hairpins("........"), [])


class WithoutViennaRNATests(unittest.TestCase):
    def test_reports_unavailable_with_reason(self):
        with mock.patch.dict(sys.modules, {"RNA": None}):
            result = rna_structure.predict_secondary_structure("GGGGAAAACCCC")
        self.assertFalse(result["available"])
        self.assertIn("ViennaRNA is not available", result["reason"])
        self.assertNotIn("structure", result)

    def test_menu_view_says_how_to_enable_it(self):
        from ui import rna_menu
        output = io.StringIO()
        with mock.patch.dict(sys.modules, {"RNA": None}), redirect_stdout(output):
            rna_menu.print_rna_secondary_structure({"rna": "GGGGAAAACCCC", "structure": {"hairpins": []}})
        self.assertIn("NOTE: ViennaRNA is not available", output.getvalue())
        self.assertIn("pip install ViennaRNA", output.getvalue())

    def test_length_limit_is_reported_not_folded(self):
        result = rna_structure.predict_secondary_structure("A" * 50, max_length=10)
        if STATUS["available"]:
            self.assertFalse(result["available"])
            self.assertIn("limited to 10 nt", result["reason"])


@NEEDS_VIENNARNA
class ViennaRNAReferenceTests(unittest.TestCase):
    def test_known_hairpin(self):
        result = rna_structure.predict_secondary_structure("GGGGAAAACCCC")
        self.assertEqual(result["structure"], "((((....))))")
        # Turner 2004, 37 °C: G-C stacks and a GAAA-closed tetraloop.
        self.assertAlmostEqual(result["mfe_kcal_per_mol"], -5.4, places=2)
        self.assertEqual(result["mfe_hairpins"][0]["stem_length"], 4)

    def test_dna_letters_and_nonstandard_bases(self):
        result = rna_structure.predict_secondary_structure("GGGGTTTTCCCCNN")
        self.assertTrue(result["available"])
        self.assertIn("2 position(s) are not A, C, G or U", result["warnings"][0])

    def test_isolated_energy_of_a_candidate(self):
        rna = "GGGGAAAACCCC"
        candidates = [c for c in find_rna_hairpins(rna) if c["stem_length"] == 4]
        result = rna_structure.predict_secondary_structure(rna, hairpin_candidates=candidates)
        item = result["hairpin_candidates"][0]
        self.assertEqual((item["pairs"], item["pairs_in_mfe"]), (4, 4))
        self.assertAlmostEqual(item["isolated_delta_g"], -5.4, places=2)

    def test_screen_against_mfe_on_a_reference_corpus(self):
        """Characterization recorded in docs/METHODS.md (RNA, hairpin screen).

        Eight 300-nt windows of pUC19 (L09137.2) and eight seeded random
        300-nt sequences. Measured with ViennaRNA 2.8.0: 64 of 74 MFE
        hairpins with stems of at least 3 bp overlap a screen candidate by
        at least half of the candidate's pairs (86 %); 1,164 of 3,307
        candidates (35 %) have a negative free energy on their own. The
        screen is sensitive but not specific; the bounds below catch a
        regression in either direction without fixing the exact counts.
        """
        lines = [line.strip() for line in PUC19.read_text().splitlines()]
        puc19 = "".join(line for line in lines if line and not line.startswith(">"))
        rng = random.Random(20261113)
        corpus = [transcribe_rna(puc19[start:start + 300]) for start in range(0, 2400, 300)]
        corpus += ["".join(rng.choice("ACGU") for _ in range(300)) for _ in range(8)]
        mfe_total = mfe_found = candidates = favorable = 0
        for rna in corpus:
            screen = find_rna_hairpins(rna)
            result = rna_structure.predict_secondary_structure(rna, hairpin_candidates=screen)
            pair_sets = []
            for candidate in screen:
                start = candidate["start"] - 1
                right = start + candidate["stem_length"] + candidate["loop_length"]
                pair_sets.append({(start + k, right + candidate["stem_length"] - 1 - k)
                                  for k in range(candidate["stem_length"])})
            for hairpin in result["mfe_hairpins"]:
                if hairpin["stem_length"] < 3:
                    continue
                inner = hairpin["start"] - 1 + hairpin["stem_length"] - 1
                stem = {(inner - k, hairpin["loop_end"] + k) for k in range(hairpin["stem_length"])}
                mfe_total += 1
                mfe_found += any(len(stem & pairs) * 2 >= len(pairs) for pairs in pair_sets)
            candidates += len(screen)
            favorable += result["candidates_favorable"]
        self.assertGreater(mfe_total, 50)
        self.assertGreaterEqual(mfe_found / mfe_total, 0.80)
        self.assertLess(favorable / candidates, 0.5)

    def test_api_entry_point(self):
        from core import api
        result = api.rna_secondary_structure("GGGGAAAACCCCAUAUGGGGAAAACCCC", sequence_type="RNA")
        self.assertEqual(result.module, "rna_structure")
        self.assertTrue(result.data["available"])
        self.assertEqual(result.parameters["tool"], "ViennaRNA")
        self.assertEqual(result.parameters["tool_version"], STATUS["version"])
        self.assertNotIn("RNA input was analyzed in the DNA alphabet (U read as T).", result.warnings)
        self.assertIn("hairpin_candidates", result.data)
        self.assertTrue(result.limitations)

    def test_menu_entry_six_draws_the_structure(self):
        import builtins
        from tests.test_terminal_ui import _scripted
        from ui import rna_menu
        buffer = io.StringIO()
        with mock.patch.object(builtins, "input", _scripted(["", "6", "", "0"])), redirect_stdout(buffer):
            rna_menu.run_rna_analysis_ui("ATGGGGAAAACCCCTAA" * 3, header="demo")
        output = buffer.getvalue()
        self.assertIn("RNA SECONDARY STRUCTURE (MFE)", output)
        self.assertIn("[6] Secondary Structure (ViennaRNA, optional)", output)
        self.assertIn("kcal/mol", output)


if __name__ == "__main__":
    unittest.main()
