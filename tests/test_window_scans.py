"""The prefix-count window scans give exactly the results of the slicing scans.

The CpG island search and the GC and complexity profiles count bases per
window. Until 2026-11-20 they sliced and recounted every window, which
took 154 s for the 4.6 Mb E. coli genome (docs/BENCHMARKS.md); they now
subtract prefix counts. The functions below are the earlier
implementations, kept verbatim in substance as the reference: on random
sequences with ambiguous and lowercase symbols, every window size and
step, the new results must be equal, float for float.
"""

import math
import random
import unittest
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import sequence as sequence_api
from analyses._sequence import cpg, profiles
from core.windows import dinucleotide_prefix_count

STANDARD = {"A", "C", "G", "T"}


def reference_cpg_windows(dna, window_size, min_gc, min_oe):
    windows = []
    for start in range(len(dna) - window_size + 1):
        window = dna[start:start + window_size]
        if any(base not in STANDARD for base in window):
            continue
        c, g, cg = window.count("C"), window.count("G"), window.count("CG")
        gc = ((c + g) / window_size) * 100
        oe = 0.0 if c == 0 or g == 0 else cg / ((c * g) / window_size)
        if gc >= min_gc and oe >= min_oe:
            windows.append((start, start + window_size - 1))
    return windows


def reference_cpg_islands(dna, window_size, min_gc, min_oe):
    """The merge and final check of find_cpg_islands, on the reference windows."""
    windows = reference_cpg_windows(dna, window_size, min_gc, min_oe)
    merged = []
    for start, end in windows:
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    islands = []
    for start, end in merged:
        region = dna[start:end + 1]
        length = len(region)
        if length < window_size:
            continue
        c, g, cg = region.count("C"), region.count("G"), region.count("CG")
        gc = ((c + g) / length) * 100
        oe = 0.0 if c == 0 or g == 0 else cg / ((c * g) / length)
        if gc < min_gc or oe < min_oe:
            continue
        islands.append({"start": start, "end": end, "start_1_based": start + 1,
                        "end_1_based": end + 1, "length": length, "gc_content": gc,
                        "cpg_count": cg, "oe_ratio": oe,
                        "criteria": {"min_length": window_size, "min_gc_content": min_gc,
                                     "min_oe_ratio": min_oe}})
    return islands


def reference_gc_profile(dna, window_size, step):
    values = []
    for start in range(0, len(dna) - window_size + 1, step):
        window = dna[start:start + window_size]
        if all(base in STANDARD for base in window):
            values.append(((window.count("G") + window.count("C")) / len(window)) * 100)
        else:
            values.append(None)
    return values


def reference_complexity(window):
    bases = [base for base in window if base in STANDARD]
    if not bases:
        return 0.0
    entropy = 0.0
    for base in sorted(STANDARD):
        count = bases.count(base)
        if count:
            probability = count / len(bases)
            entropy -= probability * math.log2(probability)
    return entropy / 2.0


def random_sequence(rng, length):
    # Mostly ACGT, CpG-rich stretches, and a few N / IUPAC / lowercase symbols.
    alphabet = "ACGT" * 20 + "CGCG" * 6 + "NRYacgt"
    return "".join(rng.choice(alphabet) for _ in range(length))


class WindowScanEquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = random.Random(20261120)
        cls.sequences = [random_sequence(rng, length) for length in (0, 1, 2, 7, 50, 199, 200, 201, 997, 3000)]

    def test_cpg_dinucleotide_prefix_matches_str_count(self):
        rng = random.Random(1)
        for dna in self.sequences:
            prefix = dinucleotide_prefix_count(dna, "CG")
            for _ in range(50):
                if len(dna) < 2:
                    break
                start = rng.randrange(len(dna))
                end = rng.randrange(start + 1, len(dna) + 1)
                self.assertEqual(prefix[end - 1] - prefix[start], dna[start:end].count("CG"))

    def test_cpg_islands_match_the_slicing_scan(self):
        for dna in self.sequences:
            for window_size, min_gc, min_oe in ((200, 50.0, 0.6), (50, 55.0, 0.65), (1, 0.0, 0.0), (7, 40.0, 0.3)):
                with self.subTest(length=len(dna), window=window_size):
                    islands = cpg.find_cpg_islands(dna, window_size, min_gc, min_oe)
                    self.assertEqual(islands, reference_cpg_islands(dna, window_size, min_gc, min_oe))

    def test_the_cpg_test_data_contains_islands(self):
        # Guard: the comparison above would be empty if no window qualified.
        self.assertTrue(any(cpg.find_cpg_islands(dna, 50, 55.0, 0.65) for dna in self.sequences))

    def test_gc_profile_matches_float_for_float(self):
        for dna in self.sequences:
            for window_size, step in ((100, 1), (100, 7), (1, 1), (13, 3)):
                with self.subTest(length=len(dna), window=window_size, step=step):
                    profile = profiles.calculate_gc_content_profile(dna, window_size, step)
                    self.assertEqual([row["gc_content"] for row in profile],
                                     reference_gc_profile(dna, window_size, step))
                    self.assertEqual([row["resolved"] for row in profile],
                                     [value is not None for value in reference_gc_profile(dna, window_size, step)])

    def test_complexity_profile_matches_float_for_float(self):
        for dna in self.sequences:
            for window_size, step in ((100, 1), (100, 9), (1, 1), (5, 2)):
                with self.subTest(length=len(dna), window=window_size, step=step):
                    profile = profiles.calculate_complexity_profile(dna, window_size, step)
                    expected = [reference_complexity(dna[s:s + window_size])
                                for s in range(0, len(dna) - window_size + 1, step)]
                    self.assertEqual([row["complexity"] for row in profile], expected)

    def test_single_window_complexity_is_unchanged(self):
        for dna in self.sequences:
            self.assertEqual(profiles.calculate_sequence_complexity(dna), reference_complexity(dna))


class LongSequenceProfileTests(unittest.TestCase):
    """Above a window limit the stored profile is thinned; the statistics are not."""

    @classmethod
    def setUpClass(cls):
        rng = random.Random(7)
        cls.sequences = [random_sequence(rng, n) for n in (99, 100, 101, 640, 2500)]

    def test_streaming_statistics_equal_the_profile_statistics(self):
        for dna in self.sequences:
            for step in (1, 3):
                with self.subTest(length=len(dna), step=step):
                    gc = profiles.calculate_gc_content_profile(dna, 100, step)
                    self.assertEqual(
                        profiles.profile_statistics_from_values(profiles.gc_window_values(dna, 100, step), 100),
                        profiles.calculate_profile_statistics(gc, "gc_content"))
                    cx = profiles.calculate_complexity_profile(dna, 100, step)
                    self.assertEqual(
                        profiles.profile_statistics_from_values(profiles.complexity_window_values(dna, 100, step), 100),
                        profiles.calculate_profile_statistics(cx, "complexity"))

    def test_thinned_profile_keeps_every_kth_window_and_exact_statistics(self):
        dna = self.sequences[-1]
        full_gc = profiles.analyze_gc_profile(dna, 100, 1)
        thin_gc = profiles.analyze_gc_profile(dna, 100, 1, max_points=500)
        stride = thin_gc["profile_stride"]
        self.assertEqual(stride, 5)                      # 2401 windows / 500, rounded up
        self.assertEqual(thin_gc["profile_windows_total"], len(full_gc["profile"]))
        self.assertEqual(thin_gc["profile"], full_gc["profile"][::stride])
        self.assertEqual(thin_gc["statistics"], full_gc["statistics"])
        full_cx = profiles.analyze_complexity_profile(dna, 100, 1)
        thin_cx = profiles.analyze_complexity_profile(dna, 100, 1, max_points=500)
        self.assertEqual(thin_cx["profile"], full_cx["profile"][::stride])
        self.assertEqual(thin_cx["statistics"], full_cx["statistics"])

    def test_below_the_limit_nothing_changes(self):
        dna = self.sequences[-1]
        self.assertEqual(profiles.analyze_gc_profile(dna, 100, 1, max_points=10_000),
                         profiles.analyze_gc_profile(dna, 100, 1))
        self.assertNotIn("profile_resolution", sequence_api.analyze_sequence(dna))

    def test_analyze_sequence_reports_the_resolution(self):
        from analyses._sequence import summary
        dna = self.sequences[-1]
        full = sequence_api.analyze_sequence(dna)
        with mock.patch.object(summary, "DEFAULT_MAX_PROFILE_POINTS", 1000):
            thinned = sequence_api.analyze_sequence(dna)
        resolution = thinned.pop("profile_resolution")
        self.assertEqual(resolution["stored_every_nth_window"], 3)
        self.assertEqual(resolution["windows_total"], 2401)
        self.assertEqual(thinned["gc_profile"], full["gc_profile"][::3])
        self.assertEqual(thinned["gc_profile_statistics"], full["gc_profile_statistics"])
        self.assertEqual(thinned["complexity_profile_statistics"], full["complexity_profile_statistics"])
        for key in set(full) - {"gc_profile", "complexity_profile"}:
            self.assertEqual(thinned[key], full[key], key)


if __name__ == "__main__":
    unittest.main()
