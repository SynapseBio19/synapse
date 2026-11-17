"""The BLAST dialog after a search, the resume path and the database manager.

Phase 7, package D. The searches are answered by a stand-in for
`NCBIBlastClient.run` that parses a recorded-shape BLAST XML document, and
the resume path by a stand-in for `NCBIBlastClient.resume`; nothing is
sent. The database-manager test builds, inspects and removes a real
`makeblastdb` database in a temporary folder and is skipped without
BLAST+.
"""

import builtins
import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from tests.test_main_session import _Script
from tests.test_p0_regressions import _blast_xml
from analyses import blast as B
from ui import blast_menu, style

DNA = "ATGGCTAGCTAGGATCCAGCTAGCTAGCTAGGCTAGCAAAGGGTTTCCCGGGAAATTTCCCGGGTTTAA"


def _fake_run(self, query, parameters):
    xml = _blast_xml(1, 32, 32, 1e-5, "A" * 32, "A" * 32)
    return B.parse_blast_xml(xml, query=query, rid="FAKE1", parameters=parameters,
                             submitted_at="2026-01-01T00:00:00+00:00")


def _fake_resume(self, rid, **_):
    query = B.BlastQuery(name="resumed_query", sequence=DNA, sequence_type="DNA")
    run = _fake_run(self, query, B.BlastParameters(program="blastn", database="nt"))
    run.metadata.warnings.append("Resumed from a stand-in in the test.")
    return run


def drive(answers, call):
    script = _Script(answers)
    output = io.StringIO()
    with mock.patch.object(builtins, "input", script), contextlib.redirect_stdout(output):
        result = call()
    return output.getvalue(), script.answers, result


# Run [1] on the selected sequence with NCBI, BLASTN, database 2 (nt), default E-value and hits.
SEARCH = ["1", "1", "2", "1", "2", "", "", "y"]


class BlastResultViewTests(unittest.TestCase):
    def test_views_before_any_search_ask_for_one(self):
        answers = ["5", "6", "7", "8", "9", "10", "11", "0"]
        output, unused, result = drive(answers, lambda: blast_menu.blast_menu(DNA, "demo"))
        self.assertEqual(unused, [])
        self.assertEqual(output.count("NOTE: Run [1] BLAST Search first."), 6)
        self.assertIn("NOTE: No BLAST hits available.", output)   # [6] Hit Detail
        self.assertIsNone(result)

    def test_every_view_after_two_searches(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            answers = (SEARCH + SEARCH
                       + ["3", "5", "6", "1", "7", "8", "9", "11", "12", "13", "14",
                          "10", str(folder / "run.json"),
                          "10", str(folder / "run.csv"),
                          "10", str(folder / "run.txt"),
                          "10", str(folder / "run.docx"),
                          "10", "",
                          "0"])
            with mock.patch.object(B.NCBIBlastClient, "run", _fake_run):
                output, unused, result = drive(answers, lambda: blast_menu.blast_menu(DNA, "demo"))
            exported = json.loads((folder / "run.json").read_text(encoding="utf-8"))
            rows = list(csv.reader((folder / "run.csv").open(encoding="utf-8")))
            text = (folder / "run.txt").read_text(encoding="utf-8")
        self.assertEqual(unused, [])
        self.assertEqual(output.count("OK: BLAST search completed."), 2)
        self.assertEqual(output.count("PRIVACY NOTICE"), 2)
        for heading in ("SEARCH CONFIGURATION", "REMOTE DATABASES", "LOCAL NCBI DATABASE CATALOG"):
            self.assertIn(heading, output)
        self.assertIn("X1", output)
        self.assertIn("ERROR: No path supplied.", output)
        self.assertIn("Use .json, .csv or .txt.", output)
        self.assertEqual(output.count("OK: Exported:"), 3)
        self.assertEqual(exported["metadata"]["rid"], "FAKE1")
        self.assertEqual(rows[1][0:1] != [], True)
        self.assertIn("X1", text)
        self.assertIsNotNone(result)

    def test_remote_search_can_be_declined(self):
        answers = ["1", "1", "2", "1", "2", "", "", "n", "0"]
        with mock.patch.object(B.NCBIBlastClient, "run", side_effect=AssertionError("must not run")):
            output, unused, result = drive(answers, lambda: blast_menu.blast_menu(DNA, "demo"))
        self.assertEqual(unused, [])
        self.assertIn("NOTE: Remote search cancelled.", output)
        self.assertIsNone(result)


class BlastResumeDialogTests(unittest.TestCase):
    def test_resume_without_rid_declined_and_accepted(self):
        answers = ["2", "",                 # no RID
                   "2", "ABCD1234", "n",    # declined
                   "2", "ABCD1234", "y",    # resumed through the stand-in
                   "5", "0"]
        with mock.patch.object(B.NCBIBlastClient, "resume", _fake_resume):
            output, unused, result = drive(answers, lambda: blast_menu.blast_menu(DNA, "demo"))
        self.assertEqual(unused, [])
        self.assertIn("ERROR: No request ID supplied.", output)
        self.assertIn("NOTE: Resume cancelled.", output)
        self.assertIn("Resumed request ABCD1234: 1 hit(s) for resumed_query", output)
        self.assertIn("NOTE: Resumed from a stand-in in the test.", output)
        self.assertIsNotNone(result)

    def test_failed_resume_is_an_error(self):
        def fail(self, rid, **_):
            raise B.BlastAnalysisError(f"NCBI does not know request {rid}.")
        with mock.patch.object(B.NCBIBlastClient, "resume", fail):
            output, unused, _ = drive(["2", "OLDRID", "y", "0"], lambda: blast_menu.blast_menu(DNA, "demo"))
        self.assertEqual(unused, [])
        self.assertIn("ERROR: NCBI does not know request OLDRID.", output)


@unittest.skipUnless(all(B.BlastExecutableManager().resolve(name) for name in ("makeblastdb", "blastdbcmd")),
                     "BLAST+ (makeblastdb, blastdbcmd) not installed")
class BlastDatabaseManagerDialogTests(unittest.TestCase):
    def test_create_inspect_verify_and_remove_a_database(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / "blastdb"
            fasta = Path(folder) / "seqs.fasta"
            fasta.write_text(">s1\n" + DNA + "\n>s2\n" + DNA[::-1] + "\n", encoding="utf-8")
            manager = B.LocalBlastDatabaseManager(root=root)
            answers = ["4", "",                                   # no FASTA path
                       "4", str(fasta), "",                       # no name
                       "4", str(fasta), "mydb", "0",              # type cancelled
                       "4", str(fasta), "mydb", "1",              # nucleotide database
                       "1",                                       # status
                       "6", "mydb",                               # information
                       "3", "mydb",                               # verify with blastdbcmd
                       "5", "absent",                             # not installed
                       "5", "mydb", "no",                         # wrong confirmation
                       "5", "mydb", "REMOVE mydb",                # removed
                       "8", "0",                                  # root change cancelled
                       "7",                                       # tools
                       "x",
                       "0"]
            output, unused, _ = drive(answers, lambda: blast_menu.blast_database_manager_menu(manager))
            remaining = list(root.rglob("mydb*"))
        self.assertEqual(unused, [])
        self.assertIn("ERROR: No FASTA path supplied.", output)
        self.assertIn("ERROR: No database name supplied.", output)
        self.assertIn("mydb", output)
        self.assertIn("ERROR: Database is not installed.", output)
        self.assertIn("NOTE: Removal cancelled.", output)
        self.assertRegex(output, r"OK: Removed \d+ files\.")
        self.assertIn("BLAST+ TOOL STATUS", output)
        self.assertIn(style.INVALID_SELECTION, output)
        # Every file makeblastdb wrote is gone, including the version 5 metadata
        # (.njs), which Remove Database left behind before Phase 7.
        self.assertEqual(remaining, [])


if __name__ == "__main__":
    unittest.main()
