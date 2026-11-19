"""The separate reverse-complement implementations agree (2026-11-19).

Seven modules have their own reverse complement (audit F6). They are kept,
because their error handling and messages differ and are part of the
menus. This test makes sure their results do not drift apart: on random
IUPAC sequences all seven equal ``core.seqtools.reverse_complement`` and,
when Biopython is installed, ``Bio.Seq.reverse_complement``. On 2026-11-19
the same comparison over 1,500 sequences found no difference.
"""

from __future__ import annotations

import random
import unittest

IUPAC = "ACGTRYSWKMBDHVN"


def _implementations():
    from analyses import cloning, mutation, primer, restriction, sequence
    from core import orf_context, seqtools
    return {
        "core.seqtools": seqtools.reverse_complement,
        "core.orf_context": orf_context.reverse_complement,
        "restriction": restriction.reverse_complement,
        "mutation (IUPAC)": mutation.reverse_complement_iupac,
        "primer": primer.reverse_complement,
        "cloning": cloning.reverse_complement,
        "sequence": sequence.get_reverse_complement,
    }


class ReverseComplementConsistencyTests(unittest.TestCase):
    def _sequences(self):
        rng = random.Random(1119)
        return ["".join(rng.choice(IUPAC) for _ in range(rng.randint(1, 60))) for _ in range(200)]

    def test_all_implementations_agree(self):
        from core.seqtools import reverse_complement
        for sequence in self._sequences():
            expected = reverse_complement(sequence)
            for name, function in _implementations().items():
                with self.subTest(implementation=name, sequence=sequence):
                    self.assertEqual(function(sequence), expected)

    def test_agreement_with_biopython(self):
        try:
            from Bio.Seq import Seq
        except ImportError:
            self.skipTest("Biopython (reference extra) is not installed")
        from core.seqtools import reverse_complement
        for sequence in self._sequences():
            self.assertEqual(reverse_complement(sequence), str(Seq(sequence).reverse_complement()))

    def test_acgt_only_variant_of_the_mutation_module(self):
        from analyses import mutation
        from core.seqtools import reverse_complement
        rng = random.Random(7)
        for _ in range(100):
            sequence = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 80)))
            self.assertEqual(mutation.reverse_complement(sequence), reverse_complement(sequence))


if __name__ == "__main__":
    unittest.main()
