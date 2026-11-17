"""The interactive session of main.py, driven in-process (Phase 7, package D).

`tests/test_main_smoke.py` runs main.py in a subprocess, which line coverage
does not follow; this file calls `main.main()` directly with scripted
answers, so the session code and the screens it draws are measured. The
"Press Enter to continue" pauses are answered automatically so that the
scripts list only real choices. The CRISPR menu entry [8] is deliberately
not driven (CRISPR work is excluded in this phase).
"""

import builtins
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from ui import style


class _Script:
    """input() replacement: pauses get Enter, every other prompt the next answer."""

    def __init__(self, answers):
        self.answers = list(answers)
        self.prompts = []

    def __call__(self, prompt=""):
        self.prompts.append(prompt)
        if style.CONTINUE_PROMPT in prompt or "Press Enter to select another sequence" in prompt:
            return ""
        if not self.answers:
            raise EOFError(f"script exhausted at prompt {prompt!r}")
        return self.answers.pop(0)


def run_session(answers):
    """Run main.main() until it exits; return (output, unused answers)."""
    import main
    script = _Script(answers)
    output = io.StringIO()
    with mock.patch.object(builtins, "input", script), contextlib.redirect_stdout(output):
        with contextlib.suppress(SystemExit):
            main.main()
    return output.getvalue(), script.answers


def write_fasta(folder, name, records):
    path = Path(folder) / name
    path.write_text("".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8")
    return str(path)


class NucleotideSessionTests(unittest.TestCase):
    def test_full_dna_session_with_reports(self):
        with tempfile.TemporaryDirectory() as folder:
            json_report = Path(folder) / "session.json"
            markdown_report = Path(folder) / "session.md"
            answers = [
                str(Path(folder) / "missing.fasta"),   # a file that does not exist
                "",                                     # then the bundled example
                "1",                                    # first record
                "1",                                    # sequence overview
                "4",                                    # protein before any ORF: refused
                "2", "5", "1",                          # ORF analysis, select ORF 1
                "3", "1", "2", "3", "4", "5", "0",      # RNA views
                "4", "1", "2", "3", "4", "5", "0",      # protein views
                "11",                                   # AI Experiment Planner: Coming Soon
                "x",                                    # invalid choice
                "P",                                    # pipeline status
                "R", str(json_report),
                "R", str(markdown_report),
                "0",
            ]
            output, unused = run_session(answers)
            report = json.loads(json_report.read_text(encoding="utf-8"))
            markdown = markdown_report.read_text(encoding="utf-8")
        self.assertEqual(unused, [])
        self.assertIn("SEQUENCE FILE ERROR", output)
        self.assertIn("No ORF has been selected.", output)
        self.assertIn("ORF CONNECTED TO PIPELINE", output)
        self.assertIn("AI Experiment Planner: COMING SOON.", output)
        self.assertIn(style.INVALID_SELECTION, output)
        self.assertRegex(output, r"ORF Analysis\s+COMPLETED \(1 stored\)")
        self.assertIn("Thank you for using Synapse.", output)
        modules = {entry["module"] for entry in report["analyses"]}
        self.assertTrue({"sequence", "orf", "rna", "protein"} <= modules, modules)
        self.assertIn("orf", markdown.lower())
        self.assertNotIn("Traceback", output)

    def test_every_module_entry_opens_and_returns(self):
        # 5 restriction/cloning, 6 primer, 7 mutation/variant, 9 BLAST, 10 literature;
        # [8] CRISPR is not driven in this phase.
        answers = ["", "1",
                   "5", "0",
                   "6", "0",
                   "7", "0",
                   "9", "0",
                   "10", "0",
                   "0"]
        output, unused = run_session(answers)
        self.assertEqual(unused, [])
        for heading in ("RESTRICTION & CLONING ANALYSIS", "PRIMER DESIGN", "MUTATION & VARIANT ANALYSIS",
                        "BLAST SEARCH", "LITERATURE SEARCH"):
            self.assertIn(heading, output)
        self.assertNotIn("ERROR:", output)

    def test_select_another_sequence_and_another_file(self):
        answers = ["", "1",
                   "S",          # back to the sequence list
                   "99", "2",    # an unlisted number, then record 2
                   "S", "F",     # back, then load another file
                   "Q"]          # quit at the file prompt
        output, unused = run_session(answers)
        self.assertEqual(unused, [])
        self.assertEqual(output.count("LOAD SEQUENCE FILE"), 2)
        self.assertIn(style.INVALID_SELECTION, output)


class RecordTypeSessionTests(unittest.TestCase):
    def test_protein_record_menu(self):
        with tempfile.TemporaryDirectory() as folder:
            fasta = write_fasta(folder, "protein.fasta", [("lysozyme_fragment", "MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV")])
            report = Path(folder) / "protein.json"
            answers = [fasta, "1",
                       "1", "1", "0",        # protein analysis: overview, back
                       "x",                  # invalid choice in the record menu
                       "R", str(report),
                       "S",                  # back to the sequence list
                       "F", "Q"]
            output, unused = run_session(answers)
            data = json.loads(report.read_text(encoding="utf-8"))
        self.assertEqual(unused, [])
        self.assertIn("PROTEIN RECORD", output)
        self.assertIn("protein (confirmed)", output)
        self.assertIn("Returning to sequence selection.", output)
        self.assertIn("protein", {entry["module"] for entry in data["analyses"]})

    def test_undetermined_type_is_declared_by_the_user(self):
        with tempfile.TemporaryDirectory() as folder:
            fasta = write_fasta(folder, "ambiguous.fasta", [("short_peptide_or_dna", "KMRSWKMRSW")])
            answers = [fasta, "1",
                       "9",          # not a type
                       "1",          # DNA: impossible for these letters? accepted as IUPAC DNA
                       "0"]
            output, unused = run_session(answers)
        self.assertEqual(unused, [])
        self.assertIn("SEQUENCE TYPE", output)
        self.assertIn("cannot be determined from its letters", output)
        self.assertIn(style.INVALID_SELECTION, output)

    def test_cancelled_type_choice_and_invalid_record(self):
        with tempfile.TemporaryDirectory() as folder:
            fasta = write_fasta(folder, "mixed.fasta", [("undetermined", "KMRSW"), ("broken", "ACGT123ACGT")])
            answers = [fasta,
                       "1", "0",     # undetermined: cancel the type question
                       "2",          # invalid letters: explained, back to the list
                       "F", "Q"]
            output, unused = run_session(answers)
        self.assertEqual(unused, [])
        self.assertIn("SEQUENCE VALIDATION ERROR", output)
        self.assertIn("Invalid character '1' at position", output)


class ReportErrorTests(unittest.TestCase):
    def test_report_error_names_the_log(self):
        import main
        output = io.StringIO()
        with mock.patch.object(builtins, "input", _Script([])), contextlib.redirect_stdout(output):
            try:
                raise ValueError("example failure")
            except ValueError as error:
                main.report_error("EXAMPLE ERROR", error)
        self.assertIn("EXAMPLE ERROR", output.getvalue())
        self.assertIn("ERROR: ValueError: example failure", output.getvalue())
        self.assertIn("Details (traceback) were written to:", output.getvalue())


if __name__ == "__main__":
    unittest.main()
