"""The faster primer screens give the results of the algorithms they replace (P2-6).

The reference functions below are the implementations as they were before
the optimization (commit bd91a24). Each optimized function must return the
same value on every input; the tiered candidate screen must return the same
candidates, in the same order, as measuring every window.

Set SYNAPSE_BENCHMARK=1 to also time a 10 kb design.
"""

from __future__ import annotations

import dataclasses
import os
import random
import time
import unittest

from analyses import primer as P
from tests._support import deterministic_dna


# ----------------------------------------------------------------------
# Reference implementations (before P2-6)
# ----------------------------------------------------------------------
def reference_reverse_complement(sequence):
    sequence = P.normalize_dna_sequence(sequence, allow_iupac=True)
    return "".join(P.IUPAC_COMPLEMENT[base] for base in reversed(sequence))


def reference_alignment_stem(left, right):
    left = P.normalize_dna_sequence(left)
    right = P.normalize_dna_sequence(right)
    right_rc = reference_reverse_complement(right)
    best = 0
    best_position = None
    for shift in range(-len(right_rc) + 1, len(left)):
        current = 0
        for index in range(len(left)):
            other_index = index - shift
            if not (0 <= other_index < len(right_rc)):
                current = 0
                continue
            if left[index] == right_rc[other_index]:
                current += 1
                if current > best:
                    best = current
                    best_position = index - current + 1
            else:
                current = 0
    return best, best_position


def reference_hairpin(sequence):
    sequence = P.normalize_dna_sequence(sequence)
    best = 0
    best_position = None
    length = len(sequence)
    max_test_stem = min(8, length // 2)
    for stem in range(max_test_stem, 2, -1):
        for left_start in range(0, length - stem - 3):
            left = sequence[left_start:left_start + stem]
            for right_start in range(left_start + stem + 3, length - stem + 1):
                right = sequence[right_start:right_start + stem]
                if left == reference_reverse_complement(right):
                    return stem, left_start
    return best, best_position


def reference_tandem_repeat(sequence, max_unit=3):
    sequence = P.normalize_dna_sequence(sequence, allow_iupac=True)
    if len(sequence) < 2:
        return 1
    best = 1
    max_unit = min(max_unit, len(sequence) // 2)
    for unit_length in range(1, max_unit + 1):
        for start in range(0, len(sequence) - 2 * unit_length + 1):
            unit = sequence[start:start + unit_length]
            copies = 1
            position = start + unit_length
            while (
                position + unit_length <= len(sequence)
                and sequence[position:position + unit_length] == unit
            ):
                copies += 1
                position += unit_length
            best = max(best, copies)
    return best


def reference_three_prime(first, second=None, window=8):
    first = P.normalize_dna_sequence(first)
    second = first if second is None else P.normalize_dna_sequence(second)
    first_terminal = first[-window:]
    second_terminal = second[-window:]
    best = 0
    for length in range(1, min(len(first_terminal), len(second_terminal)) + 1):
        if first_terminal[-length:] == reference_reverse_complement(second_terminal[:length]):
            best = length
    return best


def _random_primers(seed, count):
    rng = random.Random(seed)
    primers = []
    for _ in range(count):
        alphabet = rng.choice(["ACGT", "ACGT", "AT", "GC", "AG"])
        primers.append("".join(rng.choice(alphabet) for _ in range(rng.randint(1, 32))))
    return primers + ["GAATTC" * 4, "ACGT" * 6, "A" * 20, "GCGC" * 5, "AAGAAGAAGAAG", "acgu"]


class OptimizedHelperTests(unittest.TestCase):
    def test_reverse_complement_matches_the_lookup_on_every_iupac_symbol(self):
        rng = random.Random(3)
        for _ in range(500):
            sequence = "".join(rng.choice("ACGTRYSWKMBDHVN") for _ in range(rng.randint(1, 40)))
            self.assertEqual(P.reverse_complement(sequence), reference_reverse_complement(sequence))

    def test_dimer_stem_matches_the_diagonal_scan(self):
        primers = _random_primers(11, 600)
        for index, left in enumerate(primers):
            right = primers[(index * 7 + 3) % len(primers)]
            self.assertEqual(P._alignment_stem(left, right), reference_alignment_stem(left, right), (left, right))
            self.assertEqual(P._alignment_stem(left, left), reference_alignment_stem(left, left), left)

    def test_hairpin_matches_the_stem_by_stem_search(self):
        for primer in _random_primers(12, 1500):
            result = P.analyze_hairpin(primer)
            self.assertEqual((result.max_stem, result.best_position), reference_hairpin(primer), primer)

    def test_tandem_repeat_matches_the_copy_by_copy_count(self):
        for primer in _random_primers(13, 1500):
            for unit in (1, 2, 3, 4, 7):
                self.assertEqual(P.longest_tandem_repeat(primer, unit), reference_tandem_repeat(primer, unit), (primer, unit))

    def test_three_prime_complementarity_matches(self):
        primers = _random_primers(14, 800)
        for index, first in enumerate(primers):
            second = primers[(index * 5 + 1) % len(primers)]
            for window in (3, 8):
                self.assertEqual(P._three_prime_complementarity(first, second, window),
                                 reference_three_prime(first, second, window))
            self.assertEqual(P._three_prime_complementarity(first), reference_three_prime(first))

    def test_normalization_errors_are_unchanged_and_not_cached(self):
        for _ in range(2):
            with self.assertRaisesRegex(ValueError, "DNA sequence is empty"):
                P.normalize_dna_sequence("  ")
            with self.assertRaisesRegex(ValueError, "Invalid DNA character"):
                P.normalize_dna_sequence("ACGN")
            with self.assertRaisesRegex(TypeError, "must be a string"):
                P.normalize_dna_sequence(["A"])  # type: ignore[arg-type]
        self.assertEqual(P.normalize_dna_sequence("ACGN", allow_iupac=True), "ACGN")
        self.assertEqual(P.normalize_dna_sequence(" ac gu "), "ACGT")

    def test_duplex_cache_keeps_int_and_float_arguments_apart(self):
        as_int = P.calculate_duplex_thermodynamics("ACGTACGTACGTACGTAC", sodium_mM=50)
        as_float = P.calculate_duplex_thermodynamics("ACGTACGTACGTACGTAC", sodium_mM=50.0)
        self.assertIsInstance(as_int.sodium_mM, int)
        self.assertIsInstance(as_float.sodium_mM, float)


def _serialized(items):
    return [dataclasses.asdict(item) for item in items]


class TieredScreenTests(unittest.TestCase):
    """The tier screen keeps exactly the candidates full measurement keeps."""

    CASES = [
        ("nearest_neighbor", {}, {}),
        ("nearest_neighbor", {"min_tm": 58.0, "max_tm": 59.0, "optimal_tm": 58.5}, {"max_pair_candidates": 5}),
        ("wallace", {}, {"max_pair_candidates": 20}),
        ("salt_adjusted", {"min_tm": 55.0, "max_tm": 61.0, "optimal_tm": 58.0}, {}),
        ("santalucia", {"min_length": 15, "max_length": 24}, {"max_pair_candidates": 1}),
        ("nearest_neighbor", {"min_tm": 80.0, "max_tm": 85.0, "optimal_tm": 82.0}, {}),
    ]

    def _config(self, method, primer, extra):
        return P.ProfessionalPrimerConfig(primer=P.PrimerDesignConfig(**primer), tm_method=method, **extra)

    def test_tiered_screen_equals_measuring_every_window(self):
        rng = random.Random(5)
        for case_index, (method, primer, extra) in enumerate(self.CASES):
            sequence = "".join(rng.choice("ACGT") for _ in range(260))
            config = self._config(method, primer, extra)
            for orientation in ("FORWARD", "REVERSE"):
                start, end = (0, 260) if case_index % 2 else (20, 230)
                screened, evaluated = P.generate_professional_primer_candidates(
                    sequence, start, end, orientation, config=config)
                relaxed = P._v6_relaxed_generation_config(config.primer)
                full, full_evaluated = P._v7_measure_all_candidates(
                    sequence, start, end, orientation, config, relaxed,
                    int(config.max_pair_candidates))
                self.assertEqual(evaluated, full_evaluated)
                self.assertEqual(_serialized(screened), _serialized(full), (method, primer, extra, orientation))

    def test_upper_bound_is_never_below_the_true_score(self):
        rng = random.Random(6)
        for method in ("nearest_neighbor", "wallace", "salt_adjusted"):
            config = self._config(method, {}, {})
            P._V5_ACTIVE_CONFIG = config
            relaxed = P._v6_relaxed_generation_config(config.primer)
            sequence = "".join(rng.choice("ACGT") for _ in range(400))
            for _ in range(150):
                start = rng.randint(0, 370)
                end = start + rng.randint(18, 25)
                for orientation in ("FORWARD", "REVERSE"):
                    candidate = P._make_candidate(sequence, orientation, start, end, relaxed, config.tm_method)
                    bound = P._v7_screening_upper_bound(sequence, orientation, start, end, config, relaxed)
                    self.assertLessEqual(P._v7_screening_score(candidate, config), bound)


@unittest.skipUnless(os.environ.get("SYNAPSE_BENCHMARK"), "set SYNAPSE_BENCHMARK=1 to time a 10 kb design")
class PrimerDesignBenchmark(unittest.TestCase):
    def test_10_kb_design_time(self):
        sequence = deterministic_dna(10_000)
        started = time.process_time()
        P.design_primers_professional(sequence)
        elapsed = time.process_time() - started
        print(f"\n10 kb primer design: {elapsed:.2f} s CPU time")
        self.assertLess(elapsed, 10.0)


if __name__ == "__main__":
    unittest.main()
