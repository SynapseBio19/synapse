"""Local BLAST against a real NCBI refseq_rna database (opt-in, 2026-11-18).

Runs only with ``SYNAPSE_REFSEQ_RNA_TESTS=1`` and an installed
``refseq_rna`` (Synapse's database manager must report it as VERIFIED).
Each search over the about 84 million sequences takes about 3 minutes on
the owner's machine; the short-query case runs twice (standard search plus
rescue) and takes about 7. Nothing is downloaded or sent anywhere.

The query is nt 301-700 of NM_002046.7, Homo sapiens GAPDH transcript
variant 1, as stored in refseq_rna (read with blastdbcmd on 2026-11-18).
Several GAPDH transcript variants contain this region, so the test accepts
any human GAPDH mRNA as the top hit.

Recorded on 2026-11-18 (refseq_rna dated 30 Sep 2026):
- exact: 100 % identity and coverage, E = 0;
- three substitutions: 99.25 %;
- reverse complement: hit on the minus strand;
- 22 nt: no hit in the standard search, found by the blastn-short rescue.
"""

from __future__ import annotations

import os
from pathlib import Path
import unittest

from analyses import blast as B

GAPDH_301_700 = (
    "CATCTTCCAGGAGCGAGATCCCTCCAAAATCAAGTGGGGCGATGCTGGCGCTGAGTACGTCGTGGAGTCCACTGGCGTCTTCACCACCATGGAG"
    "AAGGCTGGGGCTCATTTGCAGGGGGGAGCCAAAAGGGTCATCATCTCTGCCCCCTCTGCTGATGCCCCCATGTTCGTCATGGGTGTGAACCATG"
    "AGAAGTATGACAACAGCCTCAAGATCATCAGCAATGCCTCCTGCACCACCAACTGCTTAGCACCCCTGGCCAAGGTCATCCATGACAACTTTGG"
    "TATCGTGGAAGGACTCATGACCACAGTCCATGCCATCACTGCCACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGC"
    "GGGGCTCTCCAGAACATCATCCCT"
)
GAPDH = "glyceraldehyde-3-phosphate dehydrogenase (GAPDH)"


def _database():
    if os.environ.get("SYNAPSE_REFSEQ_RNA_TESTS") != "1":
        return None
    database = B.LocalBlastDatabaseManager().verify_database("refseq_rna")
    return database if database.verification_status == "VERIFIED" else None


DATABASE = _database()


@unittest.skipIf(DATABASE is None, "set SYNAPSE_REFSEQ_RNA_TESTS=1 with a verified refseq_rna installed "
                                   "(about 15 minutes)")
class RefseqRnaSearchTests(unittest.TestCase):
    def _search(self, sequence, name):
        assert DATABASE is not None
        run = B.LocalBlastEngine().run(B.build_query(sequence, name=name),
                                       B.BlastParameters(program="blastn", max_target_seqs=10),
                                       DATABASE, timeout=3600)
        ranked = B.rank_hits(run.hits)
        self.assertTrue(ranked, f"no hit for {name}")
        return run, ranked[0]

    def tearDown(self):
        assert DATABASE is not None
        root = Path(DATABASE.path).resolve().parent
        self.assertEqual(sorted(p.name for p in root.rglob("synapse_blast_runtime_*")), [])

    def test_exact_fragment(self):
        _run, top = self._search(GAPDH_301_700, "gapdh_exact")
        self.assertIn(GAPDH, top.definition)
        self.assertIn("Homo sapiens", top.definition)
        self.assertEqual((top.best_identity_percent, top.best_coverage_percent), (100.0, 100.0))
        self.assertEqual(top.best_evalue, 0.0)
        self.assertNotEqual(B.calculate_evidence_score(top).evidence_band, "NOT_SIGNIFICANT")

    def test_three_substitutions(self):
        mutated = list(GAPDH_301_700)
        for position in (50, 200, 350):
            mutated[position] = {"A": "G", "G": "A", "C": "T", "T": "C"}[mutated[position]]
        _run, top = self._search("".join(mutated), "gapdh_three_substitutions")
        self.assertIn(GAPDH, top.definition)
        self.assertEqual(top.best_identity_percent, 99.25)

    def test_reverse_complement(self):
        reverse = GAPDH_301_700[::-1].translate(str.maketrans("ACGT", "TGCA"))
        _run, top = self._search(reverse, "gapdh_reverse")
        self.assertIn(GAPDH, top.definition)
        hsp = top.hsps[0]
        self.assertGreater(hsp.subject_start, hsp.subject_end)    # minus strand of the transcript

    def test_short_query_is_found_by_the_rescue(self):
        run, top = self._search(GAPDH_301_700[200:222], "gapdh_22nt")
        self.assertIn("GAPDH", top.definition)
        self.assertTrue(run.interpretation_notes[0].startswith("Standard BLASTN returned no hits"))
        self.assertEqual(top.best_identity_percent, 100.0)


if __name__ == "__main__":
    unittest.main()
