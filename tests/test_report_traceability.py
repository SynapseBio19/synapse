"""What an exported session report records (Phase 7, package H).

A session with NCBI table 11: ORF search with the bacterial code, the
selected ORF through RNA and protein analysis, JSON export. Every stored
analysis must name the software and Python version, its time, the input
and data hashes, its parameters and its limitations; RNA and protein must
state the genetic code they used, which before Phase 7 only their data
did. Nothing here is fabricated: each value is read from the result it
describes.
"""

import json
import tempfile
import unittest
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from tests.test_main_session import run_session

PUC19 = Path(__file__).parent / "fixtures" / "pUC19_L09137.2.fasta"
REQUIRED = {"module", "created_at", "synapse_version", "python_version", "input_sequence_sha256",
            "data_sha256", "parameters", "warnings", "limitations"}


class ReportTraceabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with tempfile.TemporaryDirectory() as folder:
            report = Path(folder) / "session.json"
            answers = [str(PUC19), "1",
                       "2", "6", "6", "11", "0",   # ORF settings: genetic code 11
                       "5", "1",                   # select the longest ORF
                       "3", "0",                   # RNA analysis of the ORF
                       "4", "0",                   # protein analysis of the ORF
                       "R", str(report),
                       "0"]
            cls.output, cls.unused = run_session(answers)
            cls.report = json.loads(report.read_text(encoding="utf-8"))
        cls.by_module = {entry["module"]: entry for entry in cls.report["analyses"]}

    def test_the_session_ran_as_scripted(self):
        self.assertEqual(self.unused, [])
        self.assertEqual(set(self.by_module), {"sequence", "orf", "rna", "protein"})

    def test_report_header_names_software_input_and_source_file(self):
        software, data_input = self.report["software"], self.report["input"]
        self.assertTrue({"synapse_version", "python_version", "platform"} <= set(software))
        self.assertEqual(len(data_input["sequence_sha256"]), 64)
        self.assertEqual(data_input["record"]["source_format"], "fasta")
        self.assertEqual(len(data_input["record"]["source_sha256"]), 64)
        self.assertIn("not been experimentally validated", self.report["disclaimer"])

    def test_every_analysis_carries_its_provenance(self):
        for module, entry in self.by_module.items():
            with self.subTest(module=module):
                self.assertTrue(REQUIRED <= set(entry), REQUIRED - set(entry))
                self.assertEqual(len(entry["data_sha256"]), 64)

    def test_the_genetic_code_is_recorded_wherever_it_was_used(self):
        self.assertEqual(self.by_module["orf"]["parameters"]["genetic_code_table"], 11)
        self.assertEqual(self.by_module["rna"]["parameters"]["genetic_code_table"], 11)
        self.assertEqual(self.by_module["rna"]["parameters"]["input_mode"], "orf")
        self.assertEqual(self.by_module["protein"]["parameters"]["genetic_code_table"], 11)
        self.assertIn("reading_frame", self.by_module["protein"]["parameters"])

    def test_orf_parameters_are_enough_to_repeat_the_search(self):
        parameters = self.by_module["orf"]["parameters"]
        for key in ("min_length_nt", "start_codon_policy", "start_codons", "stop_codons", "nested",
                    "ranking", "include_partial", "topology", "genetic_code_table"):
            self.assertIn(key, parameters)


if __name__ == "__main__":
    unittest.main()
