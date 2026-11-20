"""The faster RNA scans give exactly the results of the former implementations.

Changed on 2026-11-20 for genome-length input (docs/BENCHMARKS.md):

- find_rna_hairpins abandons a stem at the first mismatch beyond
  stem_length // 4, because a candidate needs paired / possible >= 0.75;
- calculate_rna_gc_profile and calculate_rna_complexity_profile take the
  window counts from prefix counts instead of slicing.

The reference functions below restate the former loops. On random RNA
with ambiguous and lowercase symbols, for several stem, loop and window
settings, the outputs must be equal, float for float and in order.
"""

import math
import random
import unittest
from typing import Any, Dict, List
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from analyses._rna import composition, structure

STANDARD = {"A", "U", "C", "G"}
WATSON_CRICK = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}


def reference_hairpins(rna, min_stem=4, min_loop=3, max_loop=12, allow_gu=False):
    structures = []
    sequence = rna.upper()
    if min_stem < 1 or min_loop < 1 or max_loop < min_loop:
        return structures
    for stem_start in range(len(sequence)):
        for stem_length in range(min_stem, min_stem + 4):
            left_end = stem_start + stem_length
            if left_end >= len(sequence):
                continue
            for loop_length in range(min_loop, max_loop + 1):
                right_start = left_end + loop_length
                right_end = right_start + stem_length
                if right_end > len(sequence):
                    continue
                left, right = sequence[stem_start:left_end], sequence[right_start:right_end]
                paired = wobble = possible = 0
                for left_base, right_base in zip(left, reversed(right)):
                    if left_base not in STANDARD or right_base not in STANDARD:
                        continue
                    possible += 1
                    if (left_base, right_base) in WATSON_CRICK:
                        paired += 1
                    elif allow_gu and {left_base, right_base} == {"G", "U"}:
                        paired += 1
                        wobble += 1
                if possible == 0:
                    continue
                fraction = paired / possible
                if fraction >= 0.75:
                    structures.append({"start": stem_start + 1, "end": right_end, "stem_length": stem_length,
                                       "loop_length": loop_length, "pairing_fraction": fraction,
                                       "wobble_pairs": wobble, "left_stem": left, "right_stem": right})
    unique, seen = [], set()
    for item in structures:
        key = (item["start"], item["end"], item["stem_length"], item["loop_length"])
        if key not in seen:
            seen.add(key)
            unique.append(item)
    unique.sort(key=lambda item: (item["pairing_fraction"], item["stem_length"]), reverse=True)
    return unique


def reference_gc(rna, window, step):
    values = []
    for start in range(0, len(rna) - window + 1, step):
        part = rna[start:start + window]
        values.append(((part.count("G") + part.count("C")) / window) * 100
                      if all(base in STANDARD for base in part) else None)
    return values


def reference_complexity(part):
    counts = [sum(1 for base in part if base == symbol) for symbol in ("A", "U", "C", "G")]
    total = sum(counts)
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts:
        if count:
            probability = count / total
            entropy -= probability * math.log2(probability)
    return entropy / 2.0


def random_rna(rng, length):
    alphabet = "ACGU" * 25 + "GCGC" * 4 + "NRYacgu"
    return "".join(rng.choice(alphabet) for _ in range(length))


class RnaScanEquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = random.Random(20261120)
        cls.sequences = [random_rna(rng, n) for n in (0, 3, 10, 25, 60, 400, 1500)]
        # A designed hairpin with one mismatch and one N in the stem.
        cls.sequences.append("GGGAACUCCAAAGGAGNUUCCCAAA")

    def test_hairpins_match_the_former_loop(self):
        for rna in self.sequences:
            settings_list: List[Dict[str, Any]] = [{}, {"allow_gu": True}, {"min_stem": 2},
                                                   {"min_stem": 6, "max_loop": 8},
                                                   {"min_stem": 9, "allow_gu": True}, {"min_loop": 0}]
            for settings in settings_list:
                with self.subTest(length=len(rna), **settings):
                    self.assertEqual(structure.find_rna_hairpins(rna, **settings),
                                     reference_hairpins(rna, **settings))

    def test_the_hairpin_data_has_candidates_with_mismatches(self):
        # Guard: the comparison must include stems that pass with a mismatch.
        found = [s for rna in self.sequences for s in reference_hairpins(rna, min_stem=6)
                 if s["pairing_fraction"] < 1.0]
        self.assertTrue(found)

    def test_gc_profile_matches_float_for_float(self):
        for rna in self.sequences:
            for window, step in ((100, 1), (100, 7), (1, 1), (9, 4)):
                with self.subTest(length=len(rna), window=window, step=step):
                    profile = composition.calculate_rna_gc_profile(rna, window, step)
                    self.assertEqual([row["gc_content"] for row in profile], reference_gc(rna, window, step))

    def test_complexity_profile_matches_float_for_float(self):
        for rna in self.sequences:
            for window, step in ((100, 1), (100, 9), (1, 1), (5, 2)):
                with self.subTest(length=len(rna), window=window, step=step):
                    profile = composition.calculate_rna_complexity_profile(rna, window, step)
                    expected = [reference_complexity(rna[s:s + window])
                                for s in range(0, len(rna) - window + 1, step)]
                    self.assertEqual([row["complexity"] for row in profile], expected)
                    self.assertEqual([composition.calculate_rna_complexity(rna[s:s + window])
                                      for s in range(0, len(rna) - window + 1, step)], expected)


class RnaLongSequenceTests(unittest.TestCase):
    def test_hairpin_count_equals_the_list_length(self):
        rng = random.Random(3)
        for rna in [random_rna(rng, n) for n in (0, 30, 800)] + ["GGGAACUCCAAAGGAGNUUCCCAAA"]:
            settings_list: List[Dict[str, Any]] = [{}, {"allow_gu": True}, {"min_stem": 6}, {"min_loop": 0}]
            for settings in settings_list:
                with self.subTest(length=len(rna), **settings):
                    self.assertEqual(structure.count_rna_hairpins(rna, **settings),
                                     len(structure.find_rna_hairpins(rna, **settings)))

    def test_thinned_rna_profiles_keep_every_kth_window(self):
        from analyses._rna import summary
        dna = random_rna(random.Random(4), 2500).replace("U", "T")
        full = summary.analyze_rna(dna)
        self.assertNotIn("profile_resolution", full)
        with mock.patch.object(summary, "DEFAULT_MAX_PROFILE_POINTS", 1000):
            thinned = summary.analyze_rna(dna)
        resolution = thinned.pop("profile_resolution")
        self.assertEqual(resolution["stored_every_nth_window"], 3)
        self.assertEqual(resolution["windows_total"], 2401)
        self.assertEqual(thinned["gc_profile"], full["gc_profile"][::3])
        self.assertEqual(thinned["complexity_profile"], full["complexity_profile"][::3])
        for key in set(full) - {"gc_profile", "complexity_profile"}:
            self.assertEqual(thinned[key], full[key], key)


if __name__ == "__main__":
    unittest.main()
