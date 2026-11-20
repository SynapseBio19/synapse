"""Per-ORF features computed faster give exactly the former values (2026-11-20).

With alternative start codons a genome-wide ORF search evaluates about
100,000 candidates per Mb, and each candidate's codon usage (twice, once
more for the diversity), ambiguity and the codon status in the scan
dominated the run (docs/BENCHMARKS.md). The functions below are the former
implementations; on random coding sequences with IUPAC symbols, every
genetic code's stop set and edge lengths, the results must be equal.
"""

import random
import unittest
from collections import Counter
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from analyses._orf import codes, features


def reference_codon_usage(coding_sequence, codon_table=None):
    stop_codons = codes._stop_codons_of(codon_table)
    codons = []
    for position in range(0, len(coding_sequence) - 2, 3):
        codon = coding_sequence[position:position + 3]
        if len(codon) < 3:
            continue
        if codon in stop_codons:
            break
        codons.append(codon)
    counts = Counter(codons)
    standard = [codon for codon in counts if all(base in "ATCG" for base in codon)]
    gc3 = [1 if codon[2] in "GC" else 0 for codon in codons if all(base in "ATCG" for base in codon)]
    return {"total_codons": len(codons), "unique_codons": len(counts), "codon_counts": dict(counts),
            "standard_codon_count": len(standard),
            "gc3_content": sum(gc3) / len(gc3) * 100 if gc3 else None}


def random_cds(rng, codons, ambiguous=False):
    alphabet = "ACGT" * 15 + ("NRYSWKM" if ambiguous else "")
    return "".join(rng.choice(alphabet) for _ in range(codons * 3 + rng.randint(0, 2)))


class OrfFeatureSpeedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = random.Random(20261120)
        cls.sequences = [""]
        for ambiguous in (False, True):
            cls.sequences += [random_cds(rng, n, ambiguous) for n in (0, 1, 2, 10, 60, 400)]
        cls.tables = [None] + [codes.codon_table_for(code) for code in (1, 2, 4, 11, 25, 27, 31)]

    def test_codon_usage_matches_the_former_loop(self):
        for sequence in self.sequences:
            for table in self.tables:
                with self.subTest(length=len(sequence)):
                    self.assertEqual(features.analyze_codon_usage(sequence, codon_table=table),
                                     reference_codon_usage(sequence, codon_table=table))

    def test_diversity_from_usage_equals_calculate_codon_diversity(self):
        for sequence in self.sequences:
            for table in self.tables:
                usage = features.analyze_codon_usage(sequence, codon_table=table)
                self.assertEqual(features.codon_diversity_from_usage(usage),
                                 features.calculate_codon_diversity(sequence, codon_table=table))

    def test_ambiguity_shortcut_equals_the_full_loop(self):
        # With the shortcut switched off the former loop runs on every input;
        # both must give the same result.
        for sequence in self.sequences:
            for table in self.tables:
                with self.subTest(length=len(sequence)):
                    fast = features.analyze_orf_ambiguity(sequence, codon_table=table)
                    with mock.patch.object(features, "_UNAMBIGUOUS", frozenset()):
                        full = features.analyze_orf_ambiguity(sequence, codon_table=table)
                    self.assertEqual(fast, full)
        self.assertTrue(any(features.analyze_orf_ambiguity(s)["ambiguous_codon_count"] for s in self.sequences))

    def test_codon_usage_stops_at_the_first_stop_of_the_table(self):
        usage = features.analyze_codon_usage("ATGTGAAAATAA", codon_table=codes.codon_table_for(1))
        self.assertEqual(usage["total_codons"], 1)
        usage = features.analyze_codon_usage("ATGTGAAAATAA", codon_table=codes.codon_table_for(4))
        self.assertEqual(usage["total_codons"], 3)          # TGA is Trp in table 4


if __name__ == "__main__":
    unittest.main()
