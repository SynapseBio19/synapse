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


class BandedAlignmentTests(unittest.TestCase):
    """Above 1,000,000 cells an exact banded alignment is tried first (P2-6)."""

    @staticmethod
    def _mutate(rng, sequence, edits):
        bases = list(sequence)
        for _ in range(edits):
            kind, position = rng.choice("sid"), rng.randrange(len(bases))
            if kind == "s":
                bases[position] = rng.choice("ACGT")
            elif kind == "i":
                bases.insert(position, rng.choice("ACGTN"))
            elif len(bases) > 1:
                del bases[position]
        return "".join(bases)

    def test_banded_and_full_alignments_are_identical(self):
        from unittest import mock
        from analyses import mutation as M
        rng = random.Random(11)
        for trial in range(150):
            length = rng.randint(20, 300)
            reference = "".join(rng.choice("ACGT") for _ in range(length))
            if trial % 10 == 0:                                   # unrelated sequences
                variant = "".join(rng.choice("ACGT") for _ in range(rng.randint(10, 300)))
            else:
                variant = self._mutate(rng, reference, rng.randint(0, length // 4 + 1))
            with self.subTest(trial=trial):
                with mock.patch.object(M, "BANDED_ALIGNMENT_MIN_CELLS", 10 ** 12):
                    full = M.global_align(reference, variant)
                with mock.patch.object(M, "BANDED_ALIGNMENT_MIN_CELLS", 0):
                    banded = M.global_align(reference, variant)
                self.assertEqual(banded, full)

    def test_the_exactness_bound(self):
        from analyses.mutation import _banded_is_exact
        # 100 x 100, half-width 32: leaving the band needs >= 66 gap columns,
        # worth at most 2 * (200 - 66) / 2 - 2 * 66 = 2 points.
        self.assertTrue(_banded_is_exact(3, 100, 100, 32, 2, -2))
        self.assertFalse(_banded_is_exact(2, 100, 100, 32, 2, -2))

    def test_banded_edit_distance_is_exact_when_within_the_band(self):
        from analyses import mutation as M
        rng = random.Random(4)
        for trial in range(120):
            reference = "".join(rng.choice("ACGTN") for _ in range(rng.randint(1, 200)))
            variant = self._mutate(rng, reference, rng.randint(0, len(reference) // 3 + 1)) or "A"
            full = M._full_edit_distance(reference, variant)
            for band in (2, 8, 32):
                band = max(band, abs(len(reference) - len(variant)))
                banded = M._banded_edit_distance(reference, variant, band)
                with self.subTest(trial=trial, band=band):
                    self.assertGreaterEqual(banded, full)          # never below the truth
                    if banded <= band:
                        self.assertEqual(banded, full)             # and exact within the band

    def test_comparing_ten_kilobases_meets_the_target(self):
        import time
        from analyses.mutation import compare_sequences
        rng = random.Random(9)
        reference = "".join(rng.choice("ACGT") for _ in range(10_000))
        variant = reference[:5000] + ("A" if reference[5000] != "A" else "C") + reference[5001:]
        started = time.perf_counter()
        result = compare_sequences(reference, variant)
        self.assertLess(time.perf_counter() - started, 10.0)     # was ~40 s (similarity matrix)
        self.assertEqual(len(result.mutations), 1)
        self.assertAlmostEqual(result.sequence_similarity_percent, 99.99)

    def test_ten_kilobases_of_similar_sequence_align_quickly(self):
        import time
        from analyses.mutation import global_align
        rng = random.Random(5)
        reference = "".join(rng.choice("ACGT") for _ in range(10_000))
        variant = self._mutate(rng, reference, 50)
        started = time.perf_counter()
        aligned_reference, aligned_variant, _ = global_align(reference, variant)
        self.assertLess(time.perf_counter() - started, 10.0)     # audit target; measured ~0.2 s
        self.assertEqual(aligned_reference.replace("-", ""), reference)
        self.assertEqual(aligned_variant.replace("-", ""), variant)


if __name__ == "__main__":
    unittest.main()
