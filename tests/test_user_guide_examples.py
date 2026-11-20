"""The statements of docs/USER_GUIDE.md about the bundled examples hold.

The guide quotes results for data/examples/pUC19.fasta (GenBank L09137.2,
the same sequence as tests/fixtures/pUC19_L09137.2.fasta with
``topology=circular`` in the header). If an analysis changes, this test
says which sentence of the guide must change with it.
"""

import unittest
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from core import api

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "data" / "examples" / "pUC19.fasta"
GUIDE = (ROOT / "docs" / "USER_GUIDE.md").read_text(encoding="utf-8")


class UserGuideExampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = api.load_records(str(EXAMPLE))[0]

    def test_example_is_pUC19_and_circular(self):
        fixture = (ROOT / "tests" / "fixtures" / "pUC19_L09137.2.fasta").read_text().splitlines()
        self.assertEqual(self.record.sequence, "".join(line.strip() for line in fixture[1:]))
        self.assertEqual(len(self.record.sequence), 2686)
        self.assertEqual(self.record.topology, "circular")
        self.assertIn("2,686 bp, circular", GUIDE)

    def test_gc_content_quoted_in_the_guide(self):
        overview = api.sequence_overview(self.record)
        self.assertEqual(f"{overview.data['gc_content']:.2f}", "50.63")
        self.assertIn("50.63 %", GUIDE)

    def test_longest_orf_is_bla(self):
        # bla (beta-lactamase) is annotated at complement(1626..2486) in L09137.2.
        orfs = api.find_orfs(self.record, genetic_code=11)
        longest = orfs.data["orfs"][0]
        self.assertEqual((longest["strand"], longest["start"] + 1, longest["stop"] + 1, longest["length"]),
                         ("-", 1626, 2486, 861))
        self.assertIn("1,626 to 2,486 (861 nt, 286 amino acids)", GUIDE)
        self.assertIn("print(longest[\"strand\"], longest[\"start\"] + 1, longest[\"stop\"] + 1)   # - 1626 2486", GUIDE)

    def test_ecori_hindiii_digest(self):
        digest = api.restriction_digest(self.record, ["EcoRI", "HindIII"])
        self.assertEqual([fragment["length"] for fragment in digest.data["fragments"]], [51, 2635])
        self.assertIn("51\nand 2,635 bp", GUIDE)
        self.assertIn("# [51, 2635]", GUIDE)

    def test_example_ships_with_the_package(self):
        self.assertIn('data = ["examples/*.fasta"]', (ROOT / "pyproject.toml").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
