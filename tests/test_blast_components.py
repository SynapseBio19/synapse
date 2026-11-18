"""Component tests for the BLAST modules (Phase 7.1).

The split of ``analyses/blast.py`` into ``analyses/_blast`` showed which
responsibilities the suite reached least: query validation (66 % of its
lines), executable discovery (62 %) and local execution (59 %). These tests
cover their documented behavior at the component boundary: every
rejected-input message, the search order for BLAST+, and the failure paths
and the short-query rescue of a local search. A fake ``subprocess.run``
stands in for BLAST+; nothing is installed or run.
"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

from analyses import blast as B
from analyses._blast import local_engine
from tests.test_p0_regressions import _blast_xml

EMPTY_XML = """<?xml version="1.0"?>
<BlastOutput><BlastOutput_program>blastn</BlastOutput_program>
<BlastOutput_iterations><Iteration><Iteration_hits></Iteration_hits></Iteration></BlastOutput_iterations>
</BlastOutput>"""


class QueryValidationTests(unittest.TestCase):
    def test_rejected_sequences_name_the_problem(self):
        cases = [
            (lambda: B.normalize_sequence(42), "Sequence must be a string."),  # type: ignore[arg-type]
            (lambda: B.normalize_sequence(" \n "), "Sequence is empty."),
            (lambda: B.normalize_sequence("ACGTXZ"), "Invalid DNA characters: X, Z"),
            (lambda: B.normalize_sequence("ACGT", "RNA"), "Invalid RNA characters: T"),
            (lambda: B.normalize_sequence("MKV1", "PROTEIN"), "Invalid protein characters: 1"),
            (lambda: B.normalize_sequence("ACGT", "PEPTIDE"), "sequence_type must be DNA, RNA or PROTEIN."),
        ]
        for call, message in cases:
            with self.subTest(message=message):
                with self.assertRaises(B.BlastAnalysisError) as raised:
                    call()
                self.assertEqual(str(raised.exception), message)

    def test_accepted_sequences_are_cleaned(self):
        self.assertEqual(B.normalize_sequence("acg u\nn", "RNA"), "ACGUN")
        self.assertEqual(B.normalize_sequence("acgt\trykm"), "ACGTRYKM")

    def test_program_and_query_type_must_match(self):
        dna = B.build_query("ACGTACGTACGT", name="dna")
        with self.assertRaises(B.BlastAnalysisError) as raised:
            B.check_program_query_type("blastp", dna)
        self.assertEqual(str(raised.exception), "blastp expects a PROTEIN query, but the query is DNA.")
        with self.assertRaises(B.BlastAnalysisError) as raised:
            B.check_program_query_type(" MegaBLAST ", dna)
        self.assertEqual(str(raised.exception), "Unsupported BLAST program: megablast")
        B.check_program_query_type(" BLASTN ", dna)
        self.assertEqual(B.programs_for_query_type("rna"), [])

    def _fasta(self, content: bytes) -> str:
        handle = tempfile.NamedTemporaryFile("wb", suffix=".fasta", delete=False)
        handle.write(content)
        handle.close()
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def test_fasta_errors(self):
        cases = [
            (b"ACGT\n>late header\nACGT\n", "Invalid FASTA file: sequence data appears before the first header."),
            (b"\n\n", "FASTA file contains no sequence records."),
            (b">one\nACGT\n>two\n", "FASTA record 2 ('two') contains no sequence."),
            (b">bad\n\xff\xfe\n", "FASTA file is not valid UTF-8:"),
        ]
        for content, message in cases:
            with self.subTest(message=message):
                with self.assertRaises(B.BlastAnalysisError) as raised:
                    B.load_fasta_queries(self._fasta(content))
                self.assertIn(message, str(raised.exception))

    def test_fasta_record_without_a_name_is_numbered(self):
        queries = B.load_fasta_queries(self._fasta(b"\xef\xbb\xbf>\nACGTACGTAC\n> second\nMKVLLAAGG\n"))
        self.assertEqual([q.name for q in queries], ["FASTA_Record_1", "second"])
        self.assertEqual([q.sequence_type for q in queries], ["DNA", "PROTEIN"])


class OrfAdapterTests(unittest.TestCase):
    def test_dicts_and_objects_are_read(self):
        self.assertEqual(B.extract_orf_context({"seq": "atg aaa", "position": "7", "id": "orf7"}),
                         ("ATGAAA", 7, "orf7"))
        self.assertEqual(B.extract_orf_context(SimpleNamespace(dna="ATGCCC", start=0, label="x")),
                         ("ATGCCC", 1, "x"))
        self.assertEqual(B.extract_orf_context({"sequence": "ATG", "start": "first"}),
                         ("ATG", 1, "Selected ORF"))

    def test_unusable_orfs_give_none(self):
        for orf in (None, {}, {"sequence": ""}, {"sequence": "ATGXQ"}):
            with self.subTest(orf=orf):
                self.assertIsNone(B.extract_orf_context(orf))

    def test_an_orf_without_a_region_is_rejected(self):
        with self.assertRaises(B.BlastAnalysisError) as raised:
            B.build_orf_query({"name": "empty"})
        self.assertIn("The selected ORF could not be extracted", str(raised.exception))


class ExecutableDiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.work = Path(tempfile.mkdtemp(prefix="synapse_bin_"))
        self.addCleanup(lambda: __import__("shutil").rmtree(self.work, ignore_errors=True))
        clean = {k: v for k, v in os.environ.items()
                 if k not in {"SYNAPSE_BLAST_BIN", "NCBI_BLAST_BIN", "BLAST_HOME", "NCBI_BLAST_HOME"}}
        patcher = mock.patch.dict(os.environ, clean, clear=True)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_search_order_explicit_then_environment(self):
        os.environ["SYNAPSE_BLAST_BIN"] = str(self.work / "synapse")
        os.environ["NCBI_BLAST_BIN"] = str(self.work / "ncbi")
        os.environ["BLAST_HOME"] = str(self.work / "home")
        bins = B.BlastExecutableManager(self.work / "explicit")._candidate_bins()
        self.assertEqual(bins[:5], [self.work / "explicit", self.work / "synapse", self.work / "ncbi",
                                    self.work / "home" / "bin", self.work / "home"])

    def test_resolve_and_require(self):
        name = "blastn.exe" if os.name == "nt" else "blastn"
        (self.work / name).write_text("")
        manager = B.BlastExecutableManager(self.work)
        with mock.patch("analyses._blast.executables.shutil.which", return_value=None):
            self.assertEqual(manager.resolve("blastn"), str(self.work / name))
            self.assertIsNone(manager.resolve("tblastx"))
            with self.assertRaises(B.BlastAnalysisError) as raised:
                manager.require("tblastx")
        self.assertEqual(str(raised.exception), "BLAST+ executable 'tblastx' was not found. "
                                                "Install NCBI BLAST+ or configure SYNAPSE_BLAST_BIN.")

    def test_update_script_next_to_a_visible_executable(self):
        (self.work / B.BLAST_UPDATE_SCRIPT).write_text("")
        blastn = self.work / "blastn"
        blastn.write_text("")

        def which(name):
            return str(blastn) if name == "blastn" else None

        with mock.patch("analyses._blast.executables.shutil.which", side_effect=which):
            found = B.BlastExecutableManager().resolve_update_script()
        self.assertEqual(Path(found or "").resolve(), (self.work / B.BLAST_UPDATE_SCRIPT).resolve())

    def test_update_script_in_a_candidate_directory(self):
        (self.work / B.BLAST_UPDATE_SCRIPT).write_text("")
        with mock.patch("analyses._blast.executables.shutil.which", return_value=None):
            self.assertEqual(B.BlastExecutableManager(self.work).resolve_update_script(),
                             str(self.work / B.BLAST_UPDATE_SCRIPT))
            self.assertIsNone(B.BlastExecutableManager(self.work / "none").resolve_update_script())


class _FakeExecutables:
    def require(self, executable):
        return executable


def _database(database_type="nucleotide"):
    return B.LocalBlastDatabase(name="mydb", database_type=database_type, path=str(Path(tempfile.gettempdir()) / "mydb"),
                                installed=True, structurally_valid=True)


class LocalEngineTests(unittest.TestCase):
    def _run(self, outputs, parameters=None, query=None, database=None, **kwargs):
        """Run the engine with a fake BLAST+ that returns ``outputs`` in turn."""
        calls = []

        def fake_run(command, **options):
            calls.append((command, options))
            result = outputs[len(calls) - 1]
            if isinstance(result, BaseException):
                raise result
            returncode, xml, stderr = result
            if xml is not None:
                Path(command[command.index("-out") + 1]).write_text(xml, encoding="utf-8")
            return subprocess.CompletedProcess(command, returncode, stdout="", stderr=stderr)

        engine = B.LocalBlastEngine(_FakeExecutables())  # type: ignore[arg-type]
        query = query or B.build_query("ACGTACGTACGTACGTACGTACGTACGTACGT", name="q")
        with mock.patch.object(local_engine.subprocess, "run", side_effect=fake_run):
            run = engine.run(query, parameters or B.BlastParameters(program="blastn"),
                             database or _database(), **kwargs)
        return run, calls

    def _error(self, *args, **kwargs):
        with self.assertRaises(B.BlastAnalysisError) as raised:
            self._run(*args, **kwargs)
        return str(raised.exception)

    def test_parameters_become_command_line_flags(self):
        parameters = B.BlastParameters(program="blastn", num_threads=2, task="megablast", word_size=11,
                                       gap_open=5, gap_extend=2, filter_query="yes")
        hit = _blast_xml(1, 32, 32, 1e-9, "A" * 32, "A" * 32)
        run, calls = self._run([(0, hit, "warning text")], parameters)
        command = calls[0][0]
        for flag, value in (("-num_threads", "2"), ("-task", "megablast"), ("-word_size", "11"),
                            ("-gapopen", "5"), ("-gapextend", "2"), ("-dust", "yes"), ("-outfmt", "5")):
            self.assertEqual(command[command.index(flag) + 1], value)
        self.assertEqual(calls[0][1]["env"]["BLASTDB"].split(os.pathsep)[0], tempfile.gettempdir())
        self.assertEqual(run.metadata.engine, "Local BLAST+")
        self.assertEqual(run.metadata.warnings[-1], "warning text")
        self.assertEqual(len(calls), 1)

    def test_protein_programs_use_seg_and_the_matrix(self):
        query = B.build_query("MKVLLAAGGHHKKLLPPQQ", name="p")
        parameters = B.BlastParameters(program="blastp", matrix="BLOSUM62", filter_query="no")
        _run, calls = self._run([(0, EMPTY_XML, "")], parameters, query, _database("protein"))
        command = calls[0][0]
        self.assertEqual(command[command.index("-seg") + 1], "no")
        self.assertEqual(command[command.index("-matrix") + 1], "BLOSUM62")

    def test_failures_are_reported(self):
        self.assertEqual(self._error([subprocess.TimeoutExpired("blastn", 5)], timeout=5),
                         "Local BLAST search exceeded the timeout of 5 seconds.")
        self.assertIn("Could not start local BLAST: denied", self._error([OSError("denied")]))
        message = self._error([(2, None, "BLAST Database error: No alias or index file found")])
        self.assertIn("Local blastn failed with exit code 2.", message)
        self.assertIn("No alias or index file found", message)
        self.assertIn("No diagnostic output was returned.", self._error([(1, None, "")]))
        self.assertEqual(self._error([(0, None, "")]),
                         "Local BLAST completed without producing an XML result file.")

    def test_preconditions_are_checked_before_running(self):
        missing = B.LocalBlastDatabase(name="gone", database_type="nucleotide", path="x")
        self.assertEqual(self._error([], database=missing),
                         "Local BLAST database 'gone' is not installed or incomplete.")
        self.assertEqual(self._error([], B.BlastParameters(program="blastp")),
                         "blastp requires a protein database, but 'mydb' is nucleotide.")
        self.assertEqual(self._error([], B.BlastParameters(program="blastn", num_threads=0)),
                         "num_threads must be at least 1.")

    def test_short_query_rescue_runs_once_when_standard_blastn_finds_nothing(self):
        hit = _blast_xml(1, 32, 32, 1e-5, "A" * 32, "A" * 32)
        run, calls = self._run([(0, EMPTY_XML, ""), (0, hit, "")])
        self.assertEqual(len(calls), 2)
        second = calls[1][0]
        self.assertEqual(second[second.index("-task") + 1], "blastn-short")
        self.assertEqual(second[second.index("-dust") + 1], "no")
        self.assertEqual(second[second.index("-evalue") + 1], calls[0][0][calls[0][0].index("-evalue") + 1])
        self.assertEqual(run.hit_count, 1)
        self.assertTrue(run.interpretation_notes[0].startswith("Standard BLASTN returned no hits"))
        self.assertTrue(run.interpretation_notes[2].startswith("Displayed hits are from the sensitive"))

    def test_rescue_without_hits_says_so(self):
        run, calls = self._run([(0, EMPTY_XML, ""), (0, EMPTY_XML, "")])
        self.assertEqual(len(calls), 2)
        self.assertEqual(run.interpretation_notes[2],
                         "No hits were found in either the standard or short-query search.")

    def test_no_rescue_with_explicit_parameters_or_long_queries(self):
        _run, calls = self._run([(0, EMPTY_XML, "")], B.BlastParameters(program="blastn", word_size=7))
        self.assertEqual(len(calls), 1)
        long_query = B.build_query("ACGT" * 60, name="long")
        _run, calls = self._run([(0, EMPTY_XML, "")], query=long_query)
        self.assertEqual(len(calls), 1)
        _run, calls = self._run([(0, EMPTY_XML, "")], allow_sensitive_retry=False)
        self.assertEqual(len(calls), 1)


class RuntimeViewCleanupTests(unittest.TestCase):
    """Synapse's temporary BLAST views do not outlive a search (2026-11-18).

    For nested multi-volume databases the local engine builds a runtime
    alias view (hard links) in ``synapse_blast_runtime_*`` next to the
    database. It was never removed. The owner's refseq_rna folder held one
    from 28 October, which the database manager then listed as a second
    81-GB database ("refseq_rna/synapse_blast_runtime_.../refseq_rna").
    """

    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="synapse_db_"))
        self.addCleanup(lambda: __import__("shutil").rmtree(self.root, ignore_errors=True))
        (self.root / "mydb").mkdir()
        self.database = B.LocalBlastDatabase(name="mydb", database_type="nucleotide",
                                             path=str(self.root / "mydb" / "mydb"),
                                             installed=True, structurally_valid=True)

    def _run(self, outcome):
        def fake_alias(database, runtime_directory):
            (Path(runtime_directory) / "mydb.nal").write_text("TITLE x\n")
            return "mydb", Path(runtime_directory)

        def fake_run(command, **options):
            if isinstance(outcome, BaseException):
                raise outcome
            Path(command[command.index("-out") + 1]).write_text(EMPTY_XML, encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        engine = B.LocalBlastEngine(_FakeExecutables())  # type: ignore[arg-type]
        query = B.build_query("ACGT" * 60, name="q")
        with mock.patch.object(local_engine, "_needs_runtime_blast_alias", return_value=True), \
                mock.patch.object(local_engine, "_create_runtime_blast_alias", side_effect=fake_alias), \
                mock.patch.object(local_engine.subprocess, "run", side_effect=fake_run):
            return engine.run(query, B.BlastParameters(program="blastn"), self.database)

    def _leftovers(self):
        return sorted(p.name for p in self.root.rglob("synapse_blast_runtime_*"))

    def test_the_view_is_removed_after_a_search(self):
        self._run(None)
        self.assertEqual(self._leftovers(), [])

    def test_the_view_is_removed_after_a_failed_search(self):
        with self.assertRaises(B.BlastAnalysisError):
            self._run(OSError("denied"))
        self.assertEqual(self._leftovers(), [])

    def test_leftover_views_are_not_listed_as_databases(self):
        # A view left behind by an interrupted search (or by Synapse before
        # this fix) must not appear as a database of its own.
        for prefix in ("synapse_blast_runtime_abc", "synapse_blast_verify_abc"):
            view = self.root / "mydb" / prefix
            view.mkdir()
            (view / "mydb.nal").write_text("TITLE view\nDBLIST mydb\n")
        (self.root / "mydb" / "mydb.nal").write_text("TITLE real\nDBLIST mydb\n")
        names = [db.name for db in B.LocalBlastDatabaseManager(self.root).list_databases()]
        self.assertIn("mydb", names)
        self.assertEqual([n for n in names if "synapse_blast_" in n], [])


class ZeroEvalueTests(unittest.TestCase):
    """BLAST reports E = 0 for its strongest matches (found 2026-11-18).

    The parser read ``Hsp_evalue`` with ``... or math.inf``, so an E-value of
    exactly 0 became infinity. A perfect 400-nt GAPDH match against the
    owner's refseq_rna was then reported as NOT_SIGNIFICANT, ranked behind
    every weaker hit and exported as the invalid JSON token ``Infinity``.
    """

    def _run(self):
        perfect = _blast_xml(1, 400, 400, "0", "A" * 400, "A" * 400)
        weaker = _blast_xml(1, 400, 400, "1e-5", "A" * 400, "A" * 400).replace("X1", "X2")
        hits = perfect.split("<Iteration_hits>")[1].split("</Iteration_hits>")[0]
        xml = weaker.replace("<Iteration_hits>", "<Iteration_hits>" + hits, 1)
        query = B.build_query("A" * 400, name="q")
        return B.parse_blast_xml(xml, query=query, rid="T", parameters=B.BlastParameters(),
                                 submitted_at="2026-01-01T00:00:00+00:00")

    def test_zero_is_kept(self):
        run = self._run()
        perfect = next(hit for hit in run.hits if hit.accession == "X1")
        self.assertEqual(perfect.best_evalue, 0.0)
        self.assertNotEqual(B.calculate_evidence_score(perfect).evidence_band, "NOT_SIGNIFICANT")

    def test_the_perfect_hit_ranks_first_and_counts_as_significant(self):
        run = self._run()
        self.assertEqual([hit.accession for hit in B.rank_hits(run.hits)], ["X1", "X2"])
        self.assertEqual(run.significant_hit_count, 2)

    def test_the_export_is_valid_json(self):
        import json
        text = B.run_to_json(self._run())
        self.assertNotIn("Infinity", text)
        json.loads(text, parse_constant=lambda token: self.fail(f"invalid JSON token {token}"))

    def test_a_missing_evalue_is_still_infinite(self):
        xml = _blast_xml(1, 50, 50, "", "A" * 50, "A" * 50)
        run = B.parse_blast_xml(xml, query=B.build_query("A" * 50, name="q"), rid="T",
                                parameters=B.BlastParameters(), submitted_at="2026-01-01T00:00:00+00:00")
        self.assertEqual(run.hits[0].best_evalue, float("inf"))


if __name__ == "__main__":
    unittest.main()
