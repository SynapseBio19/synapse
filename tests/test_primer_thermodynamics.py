"""Known-answer and cross-check tests for the nearest-neighbor Tm model.

Evidence levels:
* Primer3 manual worked example (PRIMER_TM_FORMULA=1, PRIMER_SALT_CORRECTIONS=1):
  CGTGACGTGACGGACT, 50 mM monovalent, 50 nM oligo -> dH -128.8 kcal/mol,
  salt-corrected dS -361.736 cal/(mol K), Tm 50.554 C. The manual uses
  R = 1.987; Synapse uses R = 1.98720425864083, which moves Tm by ~0.004 C.
* An independent implementation of the published equations written in this
  file (SantaLucia 1998 unified parameters typed from the paper's Table 1,
  von Ahsen 2001 divalent conversion as documented by Biopython and Primer3).
  This checks implementation conformity, not experimental accuracy.
"""

from __future__ import annotations

import math
import unittest

from tests import _support  # noqa: F401

# SantaLucia (1998) Table 1, written independently of analyses/primer.py.
UNIFIED_1998 = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
R = 1.98720425864083


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def reference_tm(seq, na=50.0, mg=1.5, dntp=0.2, oligo_nM=250.0):
    dh = ds = 0.0
    for end in (seq[0], seq[-1]):
        h, s = (0.1, -2.8) if end in "GC" else (2.3, 4.1)
        dh += h
        ds += s
    for a, b in zip(seq, seq[1:]):
        h, s = UNIFIED_1998[a + b]
        dh += h
        ds += s
    symmetric = seq == _revcomp(seq)
    if symmetric:
        ds += -1.4
    na_eq = na + (120.0 * math.sqrt(mg - dntp) if mg > dntp else 0.0)
    ds += 0.368 * (len(seq) - 1) * math.log(na_eq / 1000.0)
    c = oligo_nM * 1e-9 / (1.0 if symmetric else 4.0)
    return dh * 1000.0 / (ds + R * math.log(c)) - 273.15


PRIMERS = [
    "AGCGGATAACAATTTCACACAGGA",
    "GTAAAACGACGGCCAGT",
    "TAATACGACTCACTATAGGG",
    "CAGGAAACAGCTATGAC",
    "GCCGCGGCGCCGGCCGCG",
    "ATTATAATATTATAAT",
]


class Primer3ManualExample(unittest.TestCase):
    def test_worked_example(self):
        from analyses import primer as P
        d = P.calculate_duplex_thermodynamics("CGTGACGTGACGGACT", sodium_mM=50.0, magnesium_mM=0.0,
                                              dntp_mM=0.0, primer_nM=50.0)
        self.assertAlmostEqual(d.delta_h_kcal_mol, -128.8, places=6)
        self.assertAlmostEqual(d.delta_s_cal_mol_k, -361.736, places=3)
        self.assertAlmostEqual(d.tm_C, 50.554, delta=0.01)


class IndependentReferenceCrossCheck(unittest.TestCase):
    def test_default_model_matches_published_equations(self):
        from analyses import primer as P
        for seq in PRIMERS:
            for mg, dntp in ((0.0, 0.0), (1.5, 0.2), (3.0, 0.8), (0.5, 1.0)):
                with self.subTest(seq=seq, mg=mg, dntp=dntp):
                    observed = P.calculate_tm_nearest_neighbor(
                        seq, sodium_mM=50.0, magnesium_mM=mg, dntp_mM=dntp, primer_nM=250.0)
                    self.assertAlmostEqual(observed, reference_tm(seq, 50.0, mg, dntp, 250.0), places=6)

    def test_self_complementary_uses_total_concentration(self):
        from analyses import primer as P
        seq = "GCGCATGCGC"
        self.assertEqual(seq, _revcomp(seq))
        self.assertAlmostEqual(P.calculate_tm_nearest_neighbor(seq, magnesium_mM=0, dntp_mM=0),
                               reference_tm(seq, mg=0, dntp=0), places=6)

    def test_magnesium_raises_tm_and_dntp_chelation_lowers_it(self):
        from analyses import primer as P
        seq = PRIMERS[0]
        no_mg = P.calculate_tm_nearest_neighbor(seq, magnesium_mM=0.0, dntp_mM=0.0)
        with_mg = P.calculate_tm_nearest_neighbor(seq, magnesium_mM=1.5, dntp_mM=0.2)
        chelated = P.calculate_tm_nearest_neighbor(seq, magnesium_mM=1.5, dntp_mM=1.5)
        self.assertGreater(with_mg, no_mg + 5.0)
        self.assertAlmostEqual(chelated, no_mg, places=6)


class ParameterSetTests(unittest.TestCase):
    def test_legacy_set_reproduces_previous_values(self):
        # Values produced by Synapse before 2026-10-30 (Stage 0 golden snapshot).
        from analyses import primer as P
        expected = {"AGCGGATAACAATTTCACACAGGA": 56.519301, "GTAAAACGACGGCCAGT": 51.211379,
                    "TAATACGACTCACTATAGGG": 46.448383}
        for seq, tm in expected.items():
            observed = P.calculate_tm_nearest_neighbor(seq, nn_parameter_set="synapse_legacy",
                                                       salt_model="synapse_legacy_4x")
            self.assertAlmostEqual(observed, tm, places=5)

    def test_2004_set_differs_only_in_aa_tt_and_initiation(self):
        from analyses import primer as P
        s04 = P.NN_PARAMETER_SETS["santalucia_2004"]["stacks"]
        s98 = P.NN_PARAMETER_SETS["santalucia_1998"]["stacks"]
        self.assertEqual({k for k in s98 if s98[k] != s04[k]}, {"AA", "TT"})
        self.assertEqual(s04["AA"], (-7.6, -21.3))

    def test_unknown_names_rejected(self):
        from analyses import primer as P
        with self.assertRaises(ValueError):
            P.calculate_tm_nearest_neighbor("ACGTACGTACGT", nn_parameter_set="x")
        with self.assertRaises(ValueError):
            P.calculate_tm_nearest_neighbor("ACGTACGTACGT", salt_model="x")


class Primer3ReferenceComparison(unittest.TestCase):
    """External reference comparison with primer3-py (the ``reference`` extra).

    primer3-py wraps the Primer3 C library, an implementation that shares no
    code with Synapse. Both use SantaLucia (1998) parameters and the
    SantaLucia salt correction with the von Ahsen (2001) divalent
    conversion. Measured on 2026-11-02 with primer3-py 2.4.0: 81 primers x
    36 conditions, max |dTm| 0.0032 C; the systematic part comes from the gas
    constant (Synapse R = 1.98720 cal/(mol K), Primer3 1.987).
    """

    @unittest.skipUnless(__import__("importlib").util.find_spec("primer3"),
                         "primer3-py (reference extra) not installed")
    def test_tm_agrees_with_primer3_across_salt_conditions(self):
        import primer3
        from analyses import primer as P
        from tests._support import deterministic_dna
        template = deterministic_dna(1200)
        primers = ["CGTGACGTGACGGACT"] + [template[i:i + n] for i in range(0, 1100, 110)
                                          for n in (18, 22, 25)]
        worst = 0.0
        for seq in primers:
            for na, mg, dntp, oligo in ((50, 0, 0, 50), (50, 1.5, 0.2, 250), (100, 3.0, 0.8, 250),
                                        (50, 0.5, 0.8, 250)):
                ours = P.calculate_duplex_thermodynamics(
                    seq, sodium_mM=na, magnesium_mM=mg, dntp_mM=dntp, primer_nM=oligo).tm_C
                ref = primer3.calc_tm(seq, mv_conc=na, dv_conc=mg, dntp_conc=dntp, dna_conc=oligo,
                                      tm_method="santalucia", salt_corrections_method="santalucia")
                worst = max(worst, abs(ours - ref))
        self.assertLess(worst, 0.01)


if __name__ == "__main__":
    unittest.main()
