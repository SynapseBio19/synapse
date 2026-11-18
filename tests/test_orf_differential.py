"""Synapse's ORF finder against a naive independent finder, all NCBI codes.

Written on 2026-11-19 when all 27 NCBI genetic codes were added. The naive
finder below is deliberately simple: for every stop codon (any codon that
can terminate in the chosen code), the most upstream in-frame ATG after the
previous stop starts the ORF; both strands; complete ORFs only. With
start-codon policy "atg_only", no partial ORFs and no nested ORFs, Synapse
must report exactly the same set (strand, start, length). On 2026-11-19,
675 random sequences (25 per code) agreed; the suite runs 4 per code.
"""

from __future__ import annotations

import random
import unittest

from analyses import orf
from core.seqtools import GENETIC_CODES

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def naive_orfs(dna, stop_codons, min_length):
    found = set()
    n = len(dna)
    for strand, seq in (("+", dna), ("-", dna.translate(COMPLEMENT)[::-1])):
        for frame in range(3):
            start = None
            for i in range(frame, n - 2, 3):
                codon = seq[i:i + 3]
                if start is None and codon == "ATG":
                    start = i
                if codon in stop_codons:
                    if start is not None and i + 3 - start >= min_length:
                        found.add((strand, start if strand == "+" else n - (i + 3), i + 3 - start))
                    start = None
    return found


def synapse_orfs(dna, code, min_length):
    result = orf.analyze_orfs(dna, include_partial=False, min_length=min_length,
                              start_codon_policy="atg_only", nested=False, genetic_code=code)
    return {("+" if str(item["strand"]) in ("+", "1", "forward", "plus") else "-",
             int(item["start"]), int(item["length"])) for item in result["orfs"]}


class OrfDifferentialTests(unittest.TestCase):
    def test_every_code_agrees_with_the_naive_finder(self):
        rng = random.Random(1119)
        for code_id, code in sorted(GENETIC_CODES.items()):
            for trial in range(4):
                dna = "".join(rng.choice("ACGT") for _ in range(rng.randint(60, 900)))
                min_length = rng.choice((30, 75, 150))
                with self.subTest(code=code_id, trial=trial, length=len(dna), min_length=min_length):
                    self.assertEqual(synapse_orfs(dna, code_id, min_length),
                                     naive_orfs(dna, code.stop_codons, min_length))

    def test_the_comparison_sees_orfs(self):
        # Guard against an empty-set agreement: the fixed sequence has ORFs.
        dna = "ATG" + "GCT" * 40 + "TAA" + "CCC" * 5 + "TTA" + "AGC" * 40 + "CAT"
        found = synapse_orfs(dna, 1, 75)
        self.assertEqual(found, naive_orfs(dna, GENETIC_CODES[1].stop_codons, 75))
        self.assertEqual(len(found), 2)


if __name__ == "__main__":
    unittest.main()
