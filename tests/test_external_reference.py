"""External reference comparisons with Biopython (the ``reference`` extra).

These are external reference comparisons in the sense of the Master
Prompt (11.6): Synapse against an established implementation that shares
no code with it. They show that Synapse implements the same definitions;
they are not experimental validation. Every test is skipped, with the
reason shown, when Biopython is not installed, so the runtime stays
standard-library only. No test uses the network: Bio.Restriction carries
its REBASE-derived data locally.

Measured on 2026-11-02 with Biopython 1.88 (audit Tier 3 item 11):

* 22 enzymes: recognition site, top-strand cut and overhang length, 22/22.
* Digests of a 6,000-bp pseudo-random sequence with planted sites, linear
  and circular, three origin rotations: 132 digests, 0 differences.
* End compatibility of all 22 x 22 enzyme pairs: 484 pairs, 0 differences.
* 31 proteins: pI, GRAVY, aromaticity, extinction coefficients identical
  or within 1e-4; instability index within 2e-4 after the IG weight was
  corrected (see docs/CHANGES.md); average mass within 2e-5 relative
  (1.1e-5 measured), because the residue mass tables differ slightly
  (Synapse matches ExPASy on ubiquitin, 8564.84 Da).
"""

from __future__ import annotations

import importlib.util
import logging
import random
import unittest

from tests import _support  # noqa: F401
from tests._support import PROJECT_ROOT

HAS_BIOPYTHON = importlib.util.find_spec("Bio") is not None
SKIP_REASON = "Biopython (reference extra) not installed"


def _test_sequence():
    """6,000 bp, pseudo-random, with each rare 6-cutter planted once."""
    from data.restriction_enzymes import RESTRICTION_ENZYMES
    rng = random.Random(20261102)
    sequence = "".join(rng.choice("ACGT") for _ in range(6000))
    for name in ("EcoRI", "BamHI", "HindIII", "NdeI", "ClaI"):
        site = RESTRICTION_ENZYMES[name]["recognition_sequence"]
        position = rng.randrange(0, len(sequence) - 10)
        sequence = sequence[:position] + site + sequence[position + len(site):]
    return sequence


def _type_iis_sequence():
    """3,000 bp with every Type IIS site planted in both orientations, near both ends too."""
    from analyses.restriction import reverse_complement
    from data.restriction_enzymes import TYPE_IIS_ENZYMES
    rng = random.Random(20261111)
    sequence = "".join(rng.choice("ACGT") for _ in range(3000))
    sites = [entry["recognition_sequence"] for entry in TYPE_IIS_ENZYMES.values()]
    planted = [(site, rng.randrange(20, 2960)) for site in sites] + \
              [(reverse_complement(site), rng.randrange(20, 2960)) for site in sites] + \
              [("GGTCTC", 2990), ("GAGACC", 3)]          # cuts beyond the ends on a linear molecule
    for site, position in planted:
        sequence = sequence[:position] + site + sequence[position + len(site):]
    return sequence[:3000]


@unittest.skipUnless(HAS_BIOPYTHON, SKIP_REASON)
class TypeIISReferenceTests(unittest.TestCase):
    """Type IIS enzymes (P2-7): definitions and both-strand cuts against Biopython 1.88."""

    def test_definitions_match_rebase_data(self):
        from Bio import Restriction
        from data.restriction_enzymes import TYPE_IIS_ENZYMES
        differences = []
        for name, entry in TYPE_IIS_ENZYMES.items():
            enzyme = getattr(Restriction, name)
            ours = (entry["recognition_sequence"], entry["cut_top"], entry["cut_bottom"],
                    entry["overhang_length"])
            theirs = (str(enzyme.site), enzyme.fst5, enzyme.size + enzyme.fst3, abs(enzyme.ovhg))
            if ours != theirs:
                differences.append((name, ours, theirs))
        self.assertEqual(sorted(TYPE_IIS_ENZYMES), ["BbsI", "BsaI", "BsmBI", "SapI"])
        self.assertEqual(differences, [])

    def test_both_strand_cuts_match_on_both_topologies_and_rotated_origins(self):
        from Bio import Restriction
        from Bio.Seq import Seq
        from analyses.restriction import perform_digest
        from data.restriction_enzymes import TYPE_IIS_ENZYMES
        base = _type_iis_sequence()
        differences, cuts = [], 0
        for topology in ("linear", "circular"):
            for rotation in (0, 5, 1500, 2996):
                sequence = base[rotation:] + base[:rotation]
                for name in TYPE_IIS_ENZYMES:
                    ours = sorted(position + 1 for position in
                                  perform_digest(sequence, [name], topology=topology)["cut_positions"])
                    theirs = sorted(getattr(Restriction, name).search(
                        Seq(sequence), linear=(topology == "linear")))
                    cuts += len(theirs)
                    if ours != theirs:
                        differences.append((topology, rotation, name, ours, theirs))
        self.assertGreater(cuts, 60)
        self.assertEqual(differences, [])


@unittest.skipUnless(HAS_BIOPYTHON, SKIP_REASON)
class RestrictionReferenceTests(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.CRITICAL)

    def tearDown(self):
        logging.disable(logging.NOTSET)

    def test_enzyme_definitions_match_rebase_data(self):
        from Bio import Restriction
        from data.restriction_enzymes import RESTRICTION_ENZYMES
        differences = []
        for name, entry in RESTRICTION_ENZYMES.items():
            enzyme = getattr(Restriction, name)
            ours = (entry["recognition_sequence"], entry["cut_top"], entry["cut_bottom"],
                    entry["overhang_length"], entry["is_blunt"])
            theirs = (str(enzyme.site), enzyme.fst5, enzyme.size + enzyme.fst3,
                      abs(enzyme.ovhg), enzyme.ovhg == 0)
            if ours != theirs:
                differences.append((name, ours, theirs))
        self.assertEqual(len(RESTRICTION_ENZYMES), 22)
        self.assertEqual(differences, [])

    def test_digests_match_on_both_topologies_and_rotated_origins(self):
        from Bio import Restriction
        from Bio.Seq import Seq
        from analyses.restriction import perform_digest
        from data.restriction_enzymes import RESTRICTION_ENZYMES
        base = _test_sequence()
        differences, cuts = [], 0
        for topology in ("linear", "circular"):
            for rotation in (0, 1777, 4321):
                sequence = base[rotation:] + base[:rotation]
                for name in RESTRICTION_ENZYMES:
                    # Synapse reports the number of bases before the cut;
                    # Biopython the 1-based position of the first base after it.
                    ours = sorted(position + 1 for position in
                                  perform_digest(sequence, [name], topology=topology)["cut_positions"])
                    theirs = sorted(getattr(Restriction, name).search(
                        Seq(sequence), linear=(topology == "linear")))
                    cuts += len(theirs)
                    if ours != theirs:
                        differences.append((topology, rotation, name))
        self.assertGreater(cuts, 100)
        self.assertEqual(differences, [])

    def test_end_compatibility_matches_overhang_identity(self):
        from Bio import Restriction
        from analyses import cloning
        from data.restriction_enzymes import RESTRICTION_ENZYMES
        differences = []
        for a in RESTRICTION_ENZYMES:
            for b in RESTRICTION_ENZYMES:
                ends = [{"enzyme": name, "is_blunt": RESTRICTION_ENZYMES[name]["is_blunt"],
                         "overhang": cloning.get_overhang(name)} for name in (a, b)]
                ours = bool(cloning.ends_compatible(*ends))
                first, second = getattr(Restriction, a), getattr(Restriction, b)
                theirs = first.ovhg == second.ovhg and str(first.ovhgseq) == str(second.ovhgseq)
                if ours != theirs:
                    differences.append((a, b))
        self.assertEqual(differences, [])


@unittest.skipUnless(HAS_BIOPYTHON, SKIP_REASON)
class ProteinReferenceTests(unittest.TestCase):
    UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

    def _proteins(self):
        rng = random.Random(7)
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        return [self.UBIQUITIN] + ["M" + "".join(rng.choice(alphabet) for _ in range(rng.randrange(40, 400)))
                                   for _ in range(30)]

    def test_descriptors_match_protparam(self):
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        from analyses import protein as P
        for sequence in self._proteins():
            reference = ProteinAnalysis(sequence)
            reduced, oxidized = reference.molar_extinction_coefficient()
            extinction = P.calculate_extinction_coefficient(sequence)
            with self.subTest(length=len(sequence)):
                pairs = {
                    "pI": (P.calculate_isoelectric_point(sequence), reference.isoelectric_point(), 1e-4),
                    "GRAVY": (P.calculate_gravy(sequence), reference.gravy(), 1e-9),
                    "aromaticity": (P.calculate_aromaticity(sequence), reference.aromaticity(), 1e-9),
                    "instability": (P.calculate_instability_index(sequence),
                                    reference.instability_index(), 2e-4),
                }
                for name, (ours, theirs, tolerance) in pairs.items():
                    self.assertIsNotNone(ours, name)
                    self.assertLessEqual(abs(float(ours or 0.0) - float(theirs)), tolerance, name)
                self.assertEqual((extinction["reduced_cysteine"], extinction["oxidized_cysteine"]),
                                 (reduced, oxidized))
                mass = P.calculate_protein_molecular_weight(sequence)["known_mass"]
                self.assertLess(abs(mass - reference.molecular_weight()) / mass, 2e-5)

    def test_instability_weights_match_the_reference_table(self):
        """The dipeptide weights equal Biopython's DIWV table (ExPASy's).

        The one residual difference is FY, 33.6 here and 33.601 there.
        """
        from Bio.SeqUtils.ProtParamData import DIWV
        from analyses.protein import INSTABILITY_WEIGHTS
        differences = {first + second: (INSTABILITY_WEIGHTS[first][second], DIWV[first][second])
                       for first in DIWV for second in DIWV[first]
                       if abs(INSTABILITY_WEIGHTS[first][second] - DIWV[first][second]) > 1e-9}
        self.assertEqual(differences, {"FY": (33.6, 33.601)})


class InstabilityIgWeightTests(unittest.TestCase):
    """Known-answer test for the corrected IG weight, independent of Biopython."""

    def test_ig_dipeptide_weight_is_one(self):
        from analyses.protein import INSTABILITY_WEIGHTS, calculate_instability_index
        self.assertEqual(INSTABILITY_WEIGHTS["I"]["G"], 1.0)
        # "AIGA": dipeptides AI (1.0), IG (1.0), GA (-7.49); 10/4 * sum.
        value = calculate_instability_index("AIGA")
        self.assertIsNotNone(value)
        self.assertAlmostEqual(float(value or 0.0), 10.0 / 4 * (1.0 + 1.0 - 7.49))


# ----------------------------------------------------------------------
# pUC19 (NCBI L09137.2): the real plasmid of the 2026-10-30 audit, now a
# fixture so the comparison runs in the suite (T9, 2026-11-12).
# ----------------------------------------------------------------------
PUC19_PATH = PROJECT_ROOT / "tests" / "fixtures" / "pUC19_L09137.2.fasta"
# SHA-256 of the FASTA as NCBI E-utilities returns it (LF line ends);
# the same value was recorded at the first download on 2026-10-30.
PUC19_SHA256 = "656a4ed4597020eb58266ddd4d65f1778abff1a74655abf25ec943dcc14c152a"


def _puc19():
    text = PUC19_PATH.read_bytes().replace(b"\r\n", b"\n")
    lines = text.decode("ascii").splitlines()
    return text, lines[0], "".join(lines[1:]).upper()


class PUC19FixtureTests(unittest.TestCase):
    def test_fixture_is_the_ncbi_record(self):
        import hashlib
        text, header, sequence = _puc19()
        self.assertEqual(hashlib.sha256(text).hexdigest(), PUC19_SHA256)
        self.assertTrue(header.startswith(">L09137.2 "))
        self.assertEqual(len(sequence), 2686)
        self.assertEqual(set(sequence), set("ACGT"))


@unittest.skipUnless(HAS_BIOPYTHON, SKIP_REASON)
class PUC19RestrictionReferenceTests(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.CRITICAL)

    def tearDown(self):
        logging.disable(logging.NOTSET)

    def test_digests_of_pUC19_match_biopython(self):
        """22 enzymes, linear and circular, on pUC19 as published, rotated so an
        EcoRI site spans the origin, and rotated so a BamHI cut falls on it."""
        from Bio import Restriction
        from Bio.Seq import Seq
        from analyses.restriction import perform_digest
        from data.restriction_enzymes import RESTRICTION_ENZYMES
        _, _, puc19 = _puc19()
        ecori = puc19.find("GAATTC")
        bamhi = puc19.find("GGATCC")
        self.assertGreater(ecori, 0)
        self.assertGreater(bamhi, 0)
        rotations = (0, ecori + 3, bamhi + 1)
        differences, comparisons, cuts = [], 0, 0
        for rotation in rotations:
            sequence = puc19[rotation:] + puc19[:rotation]
            for topology in ("linear", "circular"):
                for name in RESTRICTION_ENZYMES:
                    ours = sorted(position + 1 for position in
                                  perform_digest(sequence, [name], topology=topology)["cut_positions"])
                    theirs = sorted(getattr(Restriction, name).search(
                        Seq(sequence), linear=(topology == "linear")))
                    comparisons += 1
                    cuts += len(theirs)
                    if ours != theirs:
                        differences.append((rotation, topology, name, ours, theirs))
        self.assertEqual(comparisons, 132)
        self.assertGreater(cuts, 0)
        self.assertEqual(differences, [])


if __name__ == "__main__":
    unittest.main()
