"""Primer hairpin and dimer screens against primer3-py 2.4.0 (external reference, 11.6).

Synapse classifies a structure by the length of the longest perfectly
complementary stem (and the 3' complementarity); Primer3 computes a
nearest-neighbor free energy that also weighs mismatches, loops and dangling
ends. They measure different things, so the comparison asks how far the
stem heuristic orders the thermodynamics, not whether the numbers agree.

Measured 2026-11-11 on 250 primers (150 random 18-25 nt, 50 with a planted
hairpin, 50 with a self-complementary 3' end) and 100 pairs, 50 mM Na+,
1.5 mM Mg2+, 0.2 mM dNTP, 250 nM:

* hairpin: Spearman rho(stem, dG) = -0.776; median dG LOW 0.00, MODERATE
  -1.77, HIGH -4.82 kcal/mol;
* self-dimer: rho = -0.740; median dG LOW -2.89, MODERATE -5.04, HIGH
  -11.16, VERY HIGH -9.86 kcal/mol;
* heterodimer: rho = -0.512.

The classes overlap: a MODERATE self-dimer can reach -15.6 kcal/mol, and
VERY HIGH (driven by 3' complementarity) is not more stable than HIGH in the
median. The screen is a heuristic and is documented as one (METHODS); this
test keeps the measured relationship reproducible. No network is used.
"""

from __future__ import annotations


import random
import statistics
import unittest

from tests import _support  # noqa: F401

from tests._support import primer3_unavailable

PRIMER3_UNAVAILABLE = primer3_unavailable()
CONDITIONS = dict(mv_conc=50.0, dv_conc=1.5, dntp_conc=0.2, dna_conc=250.0)


def _revcomp(sequence):
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _primer_set():
    rng = random.Random(20261111)
    primers = ["".join(rng.choice("ACGT") for _ in range(rng.randint(18, 25))) for _ in range(150)]
    for _ in range(50):
        stem = "".join(rng.choice("ACGT") for _ in range(rng.randint(4, 9)))
        loop = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 6)))
        core = stem + loop + _revcomp(stem)
        pad = 20 - len(core)
        primers.append("".join(rng.choice("ACGT") for _ in range(max(0, pad // 2))) + core
                       + "".join(rng.choice("ACGT") for _ in range(max(0, pad - pad // 2))))
    for _ in range(50):
        tail = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 6)))
        body = "".join(rng.choice("ACGT") for _ in range(20 - 2 * len(tail)))
        primers.append(body + tail + _revcomp(tail))
    return primers


def _ranks(values):
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def _spearman(x, y):
    rx, ry = _ranks(x), _ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    covariance = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    return covariance / ((sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5)


@unittest.skipIf(PRIMER3_UNAVAILABLE, PRIMER3_UNAVAILABLE or "")
class PrimerStructureReferenceTests(unittest.TestCase):
    def _measure(self, synapse, primer3_call, items):
        stems, energies, by_class = [], [], {}
        for item in items:
            ours = synapse(*item)
            theirs = primer3_call(*item, **CONDITIONS)
            energy = theirs.dg / 1000 if theirs.structure_found else 0.0
            stems.append(ours.max_stem)
            energies.append(energy)
            by_class.setdefault(ours.severity, []).append(energy)
        medians = {name: statistics.median(values) for name, values in by_class.items()}
        return _spearman(stems, energies), medians

    def test_hairpin_stems_order_the_free_energy(self):
        import primer3
        from analyses import primer
        rho, medians = self._measure(primer.analyze_hairpin, primer3.calc_hairpin,
                                     [(p,) for p in _primer_set()])
        self.assertAlmostEqual(rho, -0.776, places=3)
        self.assertGreater(medians["LOW"], medians["MODERATE"])
        self.assertGreater(medians["MODERATE"], medians["HIGH"])

    def test_self_dimer_stems_order_the_free_energy(self):
        import primer3
        from analyses import primer
        rho, medians = self._measure(primer.analyze_self_dimer, primer3.calc_homodimer,
                                     [(p,) for p in _primer_set()])
        self.assertAlmostEqual(rho, -0.740, places=3)
        self.assertGreater(medians["LOW"], medians["MODERATE"])
        self.assertGreater(medians["MODERATE"], medians["HIGH"])
        # a documented limit: 3'-driven VERY HIGH is not more stable than HIGH
        self.assertGreater(medians["VERY HIGH"], medians["HIGH"])

    def test_heterodimer_relationship_is_weaker(self):
        import primer3
        from analyses import primer
        primers = _primer_set()
        pairs = [(primers[i], primers[i + 1]) for i in range(0, 200, 2)]
        rho, _ = self._measure(primer.analyze_heterodimer, primer3.calc_heterodimer, pairs)
        self.assertAlmostEqual(rho, -0.512, places=3)


if __name__ == "__main__":
    unittest.main()
