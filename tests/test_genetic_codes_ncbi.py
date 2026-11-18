"""All NCBI genetic codes, checked against NCBI's own gc.prt (2026-11-19).

Until 2026-11-19 Synapse offered 4 of NCBI's 27 genetic codes (1, 2, 4, 11).
``core/seqtools.py`` now holds all of them, copied from gc.prt version 4.6
(``tests/fixtures/ncbi_gc/gc.prt``, see SOURCE.md there). This test reads
the file itself and compares every table: amino acids, start codons and the
codons that can terminate.

Codes 27, 28 and 31 have context-dependent codons: NCBI lists an amino acid
for them (for example TGA = W in code 27) and marks them "*" (can terminate)
in the start line. Biopython lists them only as stop codons. Synapse ends
ORFs at them and decodes them as NCBI's amino acid inside a sequence.
"""

from __future__ import annotations

import hashlib
import itertools
import re
import unittest
from pathlib import Path

from analyses import orf
from core import seqtools

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "ncbi_gc" / "gc.prt"
SHA256 = "2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04"
CODONS = ["".join(c) for c in itertools.product("TCAG", repeat=3)]


def _ncbi_tables():
    text = FIXTURE.read_text(encoding="ascii")
    tables = {}
    for entry in re.findall(r"\{(.*?)\}", text, re.S):
        ident = re.search(r"\bid\s+(\d+)", entry)
        amino = re.search(r'\bncbieaa\s+"([^"]+)"', entry)
        starts = re.search(r'\bsncbieaa\s+"([^"]+)"', entry)
        name = re.search(r'name\s+"([^"]+)"', entry)
        if ident and amino and starts and name:
            tables[int(ident.group(1))] = (" ".join(name.group(1).split()), amino.group(1), starts.group(1))
    return tables


class NcbiGeneticCodeTests(unittest.TestCase):
    def test_fixture_is_the_ncbi_file(self):
        self.assertEqual(hashlib.sha256(FIXTURE.read_bytes()).hexdigest(), SHA256)

    def test_every_ncbi_code_is_available_and_identical(self):
        ncbi = _ncbi_tables()
        self.assertEqual(sorted(ncbi), sorted(seqtools.GENETIC_CODES))
        self.assertEqual(len(ncbi), 27)
        for table_id, (name, amino, starts) in ncbi.items():
            code = seqtools.get_genetic_code(table_id)
            with self.subTest(table=table_id):
                self.assertEqual(code.name, name)
                self.assertEqual("".join(code.codons[c] for c in CODONS), amino)
                self.assertEqual(code.start_codons, {c for c, f in zip(CODONS, starts) if f == "M"})
                self.assertEqual(code.stop_codons,
                                 {c for c, a, f in zip(CODONS, amino, starts) if a == "*" or f == "*"})

    def test_context_dependent_codons(self):
        expected = {27: {"TGA"}, 28: {"TAA", "TAG", "TGA"}, 31: {"TAA", "TAG"}}
        found = {t: set(c.context_dependent_codons) for t, c in seqtools.GENETIC_CODES.items()
                 if c.context_dependent_codons}
        self.assertEqual(found, expected)
        self.assertEqual(seqtools.get_genetic_code(27).translate("ATGTGGTGAAAA", to_stop=True), "MWWK")

    def test_orfs_end_at_codons_that_can_terminate(self):
        dna = "ATG" + "AAA" * 30 + "TAA" + "GGG" * 10
        for code in (1, 28):
            with self.subTest(code=code):
                first = orf.analyze_orfs(dna, min_length=30, genetic_code=code)["orfs"][0]
                self.assertEqual((first["start"], first["stop_codon"], first["length"]), (0, "TAA", 96))

    def test_matches_biopython_except_context_dependent_codons(self):
        try:
            from Bio.Data import CodonTable
        except ImportError:
            self.skipTest("Biopython (reference extra) is not installed")
        for table_id, code in seqtools.GENETIC_CODES.items():
            reference = CodonTable.unambiguous_dna_by_id[table_id]
            with self.subTest(table=table_id):
                for codon in CODONS:
                    if codon in code.context_dependent_codons:
                        self.assertIn(codon, reference.stop_codons)
                        continue
                    expected = "*" if codon in reference.stop_codons else reference.forward_table[codon]
                    self.assertEqual(code.codons[codon], expected, codon)
                self.assertEqual(code.start_codons, set(reference.start_codons))


if __name__ == "__main__":
    unittest.main()
