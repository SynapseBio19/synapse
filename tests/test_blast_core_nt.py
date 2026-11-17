"""NCBI's core_nt database in the remote and local BLAST database lists (2026-11-17).

NCBI made core_nt (the nucleotide collection without most eukaryotic
chromosome-scale sequences) its default nucleotide database in 2024;
`nt` remains available. Synapse offers core_nt first among the remote
nucleotide databases and lists it in the update_blastdb.pl catalog.
No network: the submission is checked on the request Synapse would send.
"""

import unittest
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import blast as B


class CoreNtDatabaseTests(unittest.TestCase):
    def test_core_nt_is_a_remote_nucleotide_database_listed_first(self):
        nucleotide = [name for name, info in B.DATABASES.items() if info["type"] == "nucleotide"]
        self.assertEqual(nucleotide[0], "core_nt")
        self.assertIn("nt", nucleotide)

    def test_blastn_accepts_core_nt_and_protein_programs_refuse_it(self):
        self.assertEqual(B.BlastParameters(program="blastn", database="core_nt").normalized().database,
                         "core_nt")
        self.assertEqual(B.BlastParameters(program="tblastx", database="core_nt").normalized().database,
                         "core_nt")
        for program in ("blastp", "blastx"):
            with self.subTest(program=program), self.assertRaisesRegex(
                    B.BlastAnalysisError, "requires a protein database, but core_nt is a nucleotide"):
                B.BlastParameters(program=program, database="core_nt").normalized()

    def test_submission_names_core_nt(self):
        client = B.NCBIBlastClient()
        query = B.BlastQuery(name="q", sequence="ACGT" * 30, sequence_type="DNA")
        sent = {}

        def request(payload, method="GET"):
            sent.update(payload)
            return "RID = ABC123XYZ01R\nRTOE = 12\n"

        with mock.patch.object(client, "_request", side_effect=request):
            rid, rtoe, _ = client.submit(query, B.BlastParameters(program="blastn", database="core_nt"))
        self.assertEqual((rid, rtoe), ("ABC123XYZ01R", 12.0))
        self.assertEqual((sent["CMD"], sent["PROGRAM"], sent["DATABASE"]), ("Put", "blastn", "core_nt"))

    def test_local_catalog_lists_core_nt_without_recommending_it(self):
        entry = B.LOCAL_DATABASE_CATALOG["core_nt"]
        self.assertEqual(entry["type"], "nucleotide")
        self.assertFalse(entry["recommended"])  # very large; never downloaded automatically


if __name__ == "__main__":
    unittest.main()
