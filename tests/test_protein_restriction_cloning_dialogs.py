"""Dialog paths of the protein, restriction and cloning menus.

These three dialogs had the lowest coverage of the menus (Phase 9 report).
The scenarios drive them to the end of their branches with scripted
answers: invalid and empty input, Back at every level, each analysis and
report, and the hand-over to the pipeline. Two defects were fixed with
them (2026-11-20):

- three protein errors were printed without the ERROR prefix of the
  terminal standard;
- Back inside a cloning analysis left the cloning menu altogether instead
  of returning to it (the restriction menu returns).
"""

import random
import unittest
from pathlib import Path
from types import SimpleNamespace

from tests import _support  # noqa: F401  (blocks network and input())
from tests.test_dialog_paths import drive
from data.restriction_enzymes import RESTRICTION_ENZYMES
from ui import style

FIXTURES = Path(__file__).parent / "fixtures"


def puc19():
    lines = (FIXTURES / "pUC19_L09137.2.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines if line and not line.startswith(">"))


def enzyme_number(name):
    return str(list(RESTRICTION_ENZYMES).index(name) + 1)


def insert_sequence():
    rng = random.Random(42)
    middle = "".join(rng.choice("ACGT") for _ in range(600)).replace("GAATTC", "GAAATC").replace("AAGCTT", "AAGGTT")
    return "GAATTC" + middle + "AAGCTT"


class RecordingContext:
    def __init__(self):
        self.stored = {}

    def __getattr__(self, name):
        if name.startswith("store_"):
            def store(result, **keywords):
                self.stored[name] = (result, keywords)
            return store
        raise AttributeError(name)


def recording_pipeline():
    return SimpleNamespace(context=RecordingContext())


# ----------------------------------------------------------------------
# Protein
# ----------------------------------------------------------------------
class ProteinDialogTests(unittest.TestCase):
    def setUp(self):
        from ui import protein_menu
        self.menu = protein_menu

    def test_invalid_orf_inputs_are_errors(self):
        output, _ = drive([], lambda: self.menu.run_protein_analysis_from_orf(123))
        self.assertIn("ERROR: Protein analysis: the selected ORF sequence must be a string.", output)
        output, _ = drive([], lambda: self.menu.run_protein_analysis_from_orf("ATGXZQ"))
        self.assertIn("ERROR: Protein analysis: the selected ORF contains invalid nucleotide symbols.", output)
        output, _ = drive([], lambda: self.menu.run_protein_analysis_from_orf("TAA"))
        self.assertIn("ERROR: Protein analysis: the selected ORF produced no resolvable protein sequence.", output)

    def test_every_view_of_an_orf_and_the_pipeline_record(self):
        orf = "ATGGCTAAACGTCTGGAAGATCCGTTCTGGTAA"
        answers = []
        for choice in ("1", "2", "3", "4", "5"):
            answers += [choice, ""]
        pipeline = recording_pipeline()
        output, result = drive(answers + ["7", "0"], lambda: self.menu.run_protein_analysis_from_orf(
            orf, orf_id="ORF 1", strand="+", reading_frame="+1", orf_status="Complete", pipeline=pipeline))
        self.assertIsNone(result)
        self.assertIn(style.INVALID_SELECTION, output)
        for heading in ("PROTEIN ANALYSIS", "AMINO ACID COMPOSITION"):
            self.assertIn(heading, output.upper())
        analysis, keywords = pipeline.context.stored["store_protein_analysis"]
        self.assertEqual(keywords["protein_sequence"], "MAKRLEDPFW")
        self.assertEqual(pipeline.context.protein_source, "Selected ORF")

    def test_a_partial_orf_is_not_read_with_an_initiator_methionine(self):
        # 5'-partial: the first codon (GTG) is decoded as valine, not Met.
        pipeline = recording_pipeline()
        drive(["0"], lambda: self.menu.run_protein_analysis_from_orf(
            "GTGGCTAAATAA", orf_status="Partial - Missing Start", pipeline=pipeline))
        _, keywords = pipeline.context.stored["store_protein_analysis"]
        self.assertEqual(keywords["protein_sequence"], "VAK")
        drive(["0"], lambda: self.menu.run_protein_analysis_from_orf(
            "GTGGCTAAATAA", orf_status="Complete", pipeline=pipeline))
        _, keywords = pipeline.context.stored["store_protein_analysis"]
        self.assertEqual(keywords["protein_sequence"], "MAK")

    def test_protein_record_is_stored_and_shown(self):
        pipeline = recording_pipeline()
        output, _ = drive(["1", "", "0"], lambda: self.menu.run_protein_analysis_from_record(
            "MKT AYI\nAKQR", name="demo", pipeline=pipeline))
        self.assertEqual(pipeline.context.stored["store_protein_analysis"][1]["protein_sequence"], "MKTAYIAKQR")
        self.assertEqual(pipeline.context.protein_source, "Protein record")
        self.assertIn("demo", output)

    def test_a_record_that_is_not_protein_is_refused(self):
        with self.assertRaises(ValueError):
            drive([], lambda: self.menu.run_protein_analysis_from_record("MKT1@"))

    def test_an_empty_protein_says_so(self):
        output, _ = drive([""], lambda: self.menu.protein_analysis_menu(""))
        self.assertIn("PROTEIN ANALYSIS ERROR", output)
        self.assertIn("NOTE: No resolvable protein sequence is available.", output)


# ----------------------------------------------------------------------
# Restriction
# ----------------------------------------------------------------------
class RestrictionDialogTests(unittest.TestCase):
    def setUp(self):
        from ui import restriction_menu
        self.menu = restriction_menu
        self.dna = puc19()

    def test_every_restriction_view_with_invalid_answers(self):
        pipeline = recording_pipeline()
        answers = [
            "1",                                   # restriction analysis
            "1", "Foo", "EcoRI,HindIII", "x", "2", "",      # overview: unknown enzyme, bad topology, circular
            "2", "0",                              # map: back at the enzyme prompt
            "2", "EcoRI", "", "",                  # map: keep topology
            "3", "ALL", "", "EcoRI,BsaI", "1", "",  # digest: ALL not allowed, nothing, Type IIS
            "4", "ALL", "0",                       # report: back at the topology prompt
            "4", "ALL", "1", "",                   # report
            "7", "0",                              # invalid, back
            "x", "0",                              # invalid at the top, leave
        ]
        output, result = drive(answers, lambda: self.menu.run_restriction_cloning_ui(
            self.dna, "pUC19", sequences=[], pipeline=pipeline, topology="circular"))
        self.assertIsNone(result)
        self.assertIn("ERROR: Unknown enzyme(s): Foo", output)
        self.assertIn("ERROR: Unknown enzyme(s): ALL", output)
        self.assertIn("ERROR: Please select at least 1 enzyme(s).", output)
        self.assertGreaterEqual(output.count(style.INVALID_SELECTION), 3)
        for heading in ("RESTRICTION OVERVIEW", "RESTRICTION MAP", "DIGEST & GEL ANALYSIS",
                        "COMPLETE RESTRICTION REPORT", "TYPE IIS CUTS", "PREDICTED GEL BANDS"):
            self.assertIn(heading, output)
        stored, _ = pipeline.context.stored["store_restriction_analysis"]
        self.assertIsInstance(stored, dict)

    def test_map_without_sites_says_so(self):
        output, _ = drive(["2", "NotThere", "SpeI", "", "", "0"],
                          lambda: self.menu.run_restriction_analysis_ui("ACGT" * 50))
        self.assertIn("No restriction sites detected.", output)

    def test_an_invalid_sequence_is_reported_not_raised(self):
        output, _ = drive(["1", "EcoRI", "", "", "0"],
                          lambda: self.menu.run_restriction_analysis_ui("ACGTXX"))
        self.assertIn("ERROR: Restriction analysis failed:", output)

    def test_cloning_is_reached_and_left(self):
        # The menu does not know which record is the vector, so the loaded
        # record itself is offered as an insert; Back returns to the menu.
        sequences = [{"header": "pUC19", "sequence": self.dna}]
        output, _ = drive(["2", "0", "0"], lambda: self.menu.run_restriction_cloning_ui(
            self.dna, "pUC19", sequences=sequences))
        self.assertIn("SELECT INSERT", output)
        output, _ = drive(["2", "0"], lambda: self.menu.run_restriction_cloning_ui(
            self.dna, "pUC19", sequences=[]))
        self.assertIn("NOTE: No additional sequences are available.", output)


# ----------------------------------------------------------------------
# Cloning
# ----------------------------------------------------------------------
class CloningDialogTests(unittest.TestCase):
    def setUp(self):
        from ui import cloning_menu
        self.menu = cloning_menu
        self.vector = puc19()
        self.insert = insert_sequence()

    def run_dialog(self, answers, **keywords):
        keywords.setdefault("insert_sequence", self.insert)
        return drive(answers, lambda: self.menu.run_cloning_analysis_ui(self.vector, **keywords))

    def test_back_inside_an_analysis_returns_to_the_cloning_menu(self):
        output, result = self.run_dialog(["1", "0", "2", "0", "3", "0", "0"])
        self.assertIsNone(result)
        self.assertEqual(output.count("CLONING ANALYSIS"), 4)   # the menu is drawn again after each Back

    def test_single_enzyme_analysis(self):
        pipeline = recording_pipeline()
        output, _ = self.run_dialog(["2", "x", "1", "abc", "99", enzyme_number("EcoRI"), "", "0"],
                                    pipeline=pipeline, vector_topology="circular")
        self.assertGreaterEqual(output.count(style.INVALID_SELECTION), 3)
        self.assertIn("SINGLE-ENZYME ANALYSIS", output)
        _, keywords = pipeline.context.stored["store_cloning_analysis"]
        self.assertEqual(keywords["parameters"]["enzyme"], "EcoRI")
        self.assertEqual(keywords["parameters"]["vector_topology"], "circular")

    def test_two_enzyme_strategy_with_a_repeated_enzyme(self):
        eco, hind = enzyme_number("EcoRI"), enzyme_number("HindIII")
        pipeline = recording_pipeline()
        answers = ["3", "1",
                   eco, eco, "",                   # vector: same enzyme twice -> error, pause, again
                   eco, hind,
                   eco, hind,                      # insert
                   "x", "1", "",                   # run, report pause
                   "0"]
        output, _ = self.run_dialog(answers, pipeline=pipeline)
        self.assertIn("ERROR: The two enzymes must be different.", output)
        self.assertIn("TWO-ENZYME STRATEGY", output)
        _, keywords = pipeline.context.stored["store_cloning_analysis"]
        self.assertEqual(keywords["parameters"]["vector_enzymes"], ["EcoRI", "HindIII"])

    def test_automatic_planner_with_invalid_counts(self):
        pipeline = recording_pipeline()
        output, _ = self.run_dialog(["1", "x", "2", "abc", "0", "3", "", "0"], pipeline=pipeline)
        self.assertEqual(output.count("ERROR: Please enter a positive integer."), 2)
        self.assertIn("AUTOMATIC CLONING PLANNER", output)
        _, keywords = pipeline.context.stored["store_cloning_analysis"]
        self.assertEqual(keywords["parameters"]["top_n"], 3)

    def test_invalid_vector_and_insert(self):
        output, result = drive([], lambda: self.menu.run_cloning_analysis_ui("ATGXZ", insert_sequence="ACGT"))
        self.assertIsNone(result)
        self.assertIn("ERROR: Invalid vector sequence.", output)
        output, _ = drive([], lambda: self.menu.run_cloning_analysis_ui(self.vector, insert_sequence="ACGTJ"))
        self.assertIn("ERROR: Invalid insert sequence.", output)

    def test_insert_selection_from_the_records(self):
        sequences = [{"header": "pUC19", "sequence": self.vector},
                     {"header": "insert", "sequence": self.insert, "topology": "circular"},
                     {"header": "empty", "sequence": ""}]
        output, _ = drive(["9", "1", "0"], lambda: self.menu.run_cloning_analysis_ui(
            self.vector, sequences=sequences, vector_index=0))
        self.assertIn("SELECT INSERT", output)
        self.assertIn(style.INVALID_SELECTION, output)
        self.assertIn("(circular)", output)          # the insert keeps its record's topology
        output, _ = drive([], lambda: self.menu.run_cloning_analysis_ui(self.vector, sequences=[]))
        self.assertIn("NOTE: No additional sequences are available.", output)


if __name__ == "__main__":
    unittest.main()
