"""Phase 6 Tier 2 tests: functional defects (docs/PHASE6_INITIAL_AUDIT.md)."""

from __future__ import annotations

import unittest

from tests import _support  # noqa: F401
from tests._support import PROJECT_ROOT

LOCAL_TEST_DB = PROJECT_ROOT / "data" / "local_test" / "synapse_test"


def _blast_plus_available(*executables: str) -> bool:
    from analyses import blast as B
    manager = B.BlastExecutableManager()
    return all(manager.resolve(name) for name in executables)


# ----------------------------------------------------------------------
# D6, S7: programmatic local BLAST with a database path prefix
# ----------------------------------------------------------------------
class LocalBlastDatabasePathTests(unittest.TestCase):
    def _query(self):
        lines = (PROJECT_ROOT / "data" / "local_test" / "query.fasta").read_text().splitlines()
        return "".join(line.strip() for line in lines if not line.startswith(">"))

    @unittest.skipUnless(_blast_plus_available("blastn", "blastdbcmd"),
                         "BLAST+ (blastn, blastdbcmd) not installed")
    def test_prefix_path_runs_and_records_provenance(self):
        """External tool execution: the bundled database addressed by its prefix."""
        from analyses import blast as B
        run = B.run_blast_analysis(self._query(), engine="local", program="blastn",
                                   local_database_path=str(LOCAL_TEST_DB))
        self.assertEqual([(hit.hit_id, round(hit.best_coverage_percent, 2)) for hit in run.hits],
                         [("Test_Gene_001", 100.0), ("Test_Gene_002", 98.04)])  # audit section 6
        details = run.metadata.database_details
        self.assertEqual(details["title"], "Synapse Local Test Database")
        self.assertEqual(details["sequence_count"], 3)
        self.assertIsNotNone(details["database_date"])

    def test_missing_prefix_is_reported(self):
        from analyses import blast as B
        with self.assertRaisesRegex(B.BlastAnalysisError, "No BLAST database found at prefix"):
            B.run_blast_analysis(self._query(), engine="local", program="blastn",
                                 local_database_path=str(LOCAL_TEST_DB.parent / "missing_db"))

    @unittest.skipUnless(_blast_plus_available("blastdbcmd"), "BLAST+ (blastdbcmd) not installed")
    def test_program_database_type_mismatch_is_refused(self):
        from analyses import blast as B
        with self.assertRaisesRegex(B.BlastAnalysisError, "requires a protein database"):
            B.run_blast_analysis("MKTAYIAKQRQISFVKSHFSRQ", engine="local", program="blastp",
                                 local_database_path=str(LOCAL_TEST_DB))


if __name__ == "__main__":
    unittest.main()
