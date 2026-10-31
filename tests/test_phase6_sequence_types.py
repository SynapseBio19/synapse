"""Phase 6 tests: sequence types across modules (audit D14, D13)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tests import _support  # noqa: F401


PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ"


def _blast_plus_available(*executables: str) -> bool:
    from analyses import blast as B
    manager = B.BlastExecutableManager()
    return all(manager.resolve(name) for name in executables)


# ----------------------------------------------------------------------
# D14-F2: one query-type vocabulary in BLAST
# ----------------------------------------------------------------------
class BlastQueryTypeVocabularyTests(unittest.TestCase):
    def test_program_table_uses_query_vocabulary(self):
        from analyses import blast as B
        for name, info in B.PROGRAMS.items():
            with self.subTest(program=name):
                self.assertIn(info["query_type"], B.QUERY_TYPES)
        self.assertEqual(B.build_query(PROTEIN).sequence_type, "PROTEIN")
        self.assertEqual(B.build_query("ACGTACGT").sequence_type, "DNA")

    def test_programs_offered_for_each_query_type(self):
        from analyses import blast as B
        self.assertEqual(B.programs_for_query_type("PROTEIN"), ["blastp", "tblastn"])
        self.assertEqual(B.programs_for_query_type("DNA"), ["blastn", "blastx", "tblastx"])

    def test_check_accepts_matching_and_rejects_other_type(self):
        from analyses import blast as B
        protein, dna = B.build_query(PROTEIN), B.build_query("ACGTACGTACGT")
        for program in ("blastp", "tblastn"):
            B.check_program_query_type(program, protein)
            with self.assertRaisesRegex(B.BlastAnalysisError, "expects a PROTEIN query"):
                B.check_program_query_type(program, dna)
        for program in ("blastn", "blastx", "tblastx"):
            B.check_program_query_type(program, dna)
            with self.assertRaisesRegex(B.BlastAnalysisError, "expects a DNA query"):
                B.check_program_query_type(program, protein)

    def test_remote_run_rejects_wrong_query_type_before_submitting(self):
        from analyses import blast as B
        client = B.NCBIBlastClient()
        submitted = []
        client.submit = lambda *args, **kwargs: submitted.append(args)  # type: ignore[method-assign]
        with self.assertRaisesRegex(B.BlastAnalysisError, "expects a PROTEIN query"):
            client.run(B.build_query("ACGTACGTACGT"), B.BlastParameters(program="blastp", database="nr"))
        self.assertEqual(submitted, [])

    @unittest.skipUnless(_blast_plus_available("makeblastdb", "blastp"),
                         "BLAST+ (makeblastdb, blastp) not installed")
    def test_local_blastp_runs_protein_query(self):
        """External tool execution: a protein query against a two-record protein database."""
        from analyses import blast as B
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            fasta = Path(tmp) / "proteins.fasta"
            fasta.write_text(f">target1\n{PROTEIN}\n>target2\n{'GSGSWWPPHHCC' * 6}\n")
            database = B.LocalBlastDatabaseManager(root=tmp).create_custom_database(
                fasta, "mini_protein", "protein")
            run = B.LocalBlastEngine().run(
                B.build_query(PROTEIN[5:45], name="Q"),
                B.BlastParameters(program="blastp", database="nr"),
                database,
            )
        self.assertEqual([hit.hit_id for hit in run.hits][:1], ["target1"])


if __name__ == "__main__":
    unittest.main()
