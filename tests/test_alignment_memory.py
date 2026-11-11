"""P2-6: the mutation alignment keeps one traceback byte per cell.

The global alignment held two (n+1) x (m+1) matrices of Python objects,
measured at about 46 bytes per cell (184 MB at 2,000 nt, ~4.5 GB at
10,000 nt). It now keeps two score rows and a bytearray traceback with the
same recurrence and tie-breaking; 49 alignments (lengths 1-400, indels, IUPAC
symbols, refused empty input) were byte-identical before and after.
"""

from __future__ import annotations

import random
import tracemalloc
import unittest

from tests import _support  # noqa: F401


class AlignmentMemoryTests(unittest.TestCase):
    def test_peak_memory_is_about_one_byte_per_cell(self):
        from analyses.mutation import global_align
        rng = random.Random(3)
        reference = "".join(rng.choice("ACGT") for _ in range(600))
        variant = reference[:200] + "GG" + reference[203:]
        tracemalloc.start()
        try:
            global_align(reference, variant)
            peak = tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()
        cells = (len(reference) + 1) * (len(variant) + 1)
        self.assertLess(peak / cells, 3.0)          # was ~46

    def test_tie_breaking_prefers_diagonal_then_gap_in_the_variant(self):
        from analyses.mutation import global_align
        self.assertEqual(global_align("AAAA", "AA"), ("AAAA", "--AA", 0))
        self.assertEqual(global_align("ACACAC", "ACAC"), ("ACACAC", "--ACAC", 4))
        self.assertEqual(global_align("RYKM", "ACGT"), ("RYKM", "ACGT", 5))

    def test_empty_sequences_are_refused_as_before(self):
        from analyses.mutation import MutationAnalysisError, global_align
        for reference, variant in (("", "ACG"), ("A", "")):
            with self.assertRaises(MutationAnalysisError):
                global_align(reference, variant)

    def test_iupac_symbols_match_their_bases(self):
        from analyses.mutation import global_align
        self.assertEqual(global_align("ACGTNACGT", "ACGTAACGT"), ("ACGTNACGT", "ACGTAACGT", 18))


if __name__ == "__main__":
    unittest.main()
