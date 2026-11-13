"""ORF search compared with NCBI ORFfinder (web service, 2026-11-13).

The reference is ORFfinder's result table for pUC19 (L09137.2, entered as
an accession) with genetic code 11, four settings:
tests/fixtures/orffinder/pUC19_L09137.2_orffinder.json. ORFfinder was run
live once; the table is stored so the comparison is reproducible offline.

Finding: every ORFfinder ORF is reported by Synapse with the same strand,
start, stop, length and partial status. The only difference is the length
threshold convention: ORFfinder's "minimal ORF length" does not count the
stop codon, Synapse's ``min_length`` does. ``min_length = L + 3`` in
Synapse reproduces ORFfinder's "L" exactly.
"""

import hashlib
import json
import unittest
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from core import api

FIXTURES = Path(__file__).parent / "fixtures"
REFERENCE = FIXTURES / "orffinder" / "pUC19_L09137.2_orffinder.json"
REFERENCE_SHA256 = "2053f818059ad16849b2444aceaa7b9f21cf18b0bacb6337f6176cdacf599e11"


def puc19():
    lines = (FIXTURES / "pUC19_L09137.2.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines if line and not line.startswith(">"))


def orffinder_keys(run):
    keys = set()
    for orf in run["orfs"]:
        start, stop = orf["start"], orf["stop"]
        partial = stop.startswith(">")
        stop_value = int(stop.lstrip(">"))
        low, high = (int(start), stop_value) if orf["strand"] == "+" else (stop_value, int(start))
        keys.add((orf["strand"], low, high, orf["length_nt"], partial))
    return keys


def synapse_keys(run, min_length):
    result = api.find_orfs(puc19(), min_length=min_length, start_codon_policy=run["synapse_policy"],
                           genetic_code=11, include_partial=True, top_n=10_000)
    # Synapse records: 0-based start and 0-based last base; length includes the stop codon.
    return {(o["strand"], o["start"] + 1, o["stop"] + 1, o["length"], o["status"] != "Complete")
            for o in result.data["orfs"]}


class OrffinderComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = json.loads(REFERENCE.read_text())

    def test_reference_file_is_unchanged(self):
        digest = hashlib.sha256(REFERENCE.read_bytes()).hexdigest()
        self.assertEqual(digest, REFERENCE_SHA256)
        self.assertEqual(sum(len(run["orfs"]) for run in self.reference["runs"].values()), 69)

    def test_same_orfs_when_the_threshold_excludes_the_stop_codon(self):
        for name, run in self.reference["runs"].items():
            with self.subTest(run=name):
                self.assertEqual(synapse_keys(run, run["minimal_length_nt"] + 3), orffinder_keys(run))

    def test_with_the_same_number_only_orfs_of_exactly_that_length_differ(self):
        for name, run in self.reference["runs"].items():
            with self.subTest(run=name):
                synapse = synapse_keys(run, run["minimal_length_nt"])
                reference = orffinder_keys(run)
                self.assertLessEqual(reference, synapse)
                self.assertTrue(all(key[3] == run["minimal_length_nt"] for key in synapse - reference))

    def test_alternative_starts_and_partial_orfs_are_covered(self):
        run = self.reference["runs"]["alternative_min75"]
        self.assertEqual(len(run["orfs"]), 44)
        self.assertEqual(sum(1 for key in orffinder_keys(run) if key[4]), 2)
        result = api.find_orfs(puc19(), min_length=78, start_codon_policy="atg_and_alternative",
                               genetic_code=11, include_partial=True, top_n=10_000)
        starts = {o["start_codon"] for o in result.data["orfs"]}
        self.assertEqual(starts, {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"})


if __name__ == "__main__":
    unittest.main()
