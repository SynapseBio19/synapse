"""HGVS DNA descriptions compared with Mutalyzer 3 (LUMC), recorded 2026-11-13.

Mutalyzer's ``description_extract`` endpoint (API 2.2, Mutalyzer
3.2.0.dev0, https://mutalyzer.nl/api/description_extract/) computes the
canonical HGVS description from a reference and an observed sequence,
applying the 3' rule and recognizing duplications. It was run once with
the synthetic sequences below; its answers are recorded here so the test
runs offline. Mutalyzer writes a duplicated mononucleotide pair as a
repeat (``9_10T[4]``); ``9_10dup`` describes the same sequence and is
equally valid HGVS.

The protein forms Synapse writes (``p.(R2C)``, ``p.(W24*)``, ``p.(M1?)``,
``p.(R12=)``, ``p.(*110Qext*?)``) were parsed by Mutalyzer's
``description_to_model`` without error and flagged as predicted; a
truncated form (``p.(R2C``) was rejected (ESYNTAXUEOF). Recorded in
docs/METHODS.md.
"""

import unittest

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import mutation as M

REFERENCE = "ATGAAACGTTGGCTGCAAGGTTAA"
# name: (observed sequence, Mutalyzer description, Synapse description if it differs in form only)
MUTALYZER = {
    "substitution": ("ATGAAATGTTGGCTGCAAGGTTAA", "7C>T", None),
    "deletion_in_AAA_run": ("ATGAACGTTGGCTGCAAGGTTAA", "6del", None),
    "deletion_single": ("ATGAAACGTTGGCTGAAGGTTAA", "16del", None),
    "deletion_three": ("ATGAAACGTCTGCAAGGTTAA", "10_12del", None),
    "duplication_TT": ("ATGAAACGTTTTGGCTGCAAGGTTAA", "9_10T[4]", "9_10dup"),
    "duplication_A": ("ATGAAAACGTTGGCTGCAAGGTTAA", "6dup", None),
    "duplication_after_shift": ("ATGAAACGTTGGCGCTGCAAGGTTAA", "12_13dup", None),
    "duplication_G": ("ATGAAACGTTGGGCTGCAAGGTTAA", "12dup", None),
    "adjacent_substitutions": ("ATGAAACGTTGGCTGCAAGCCTAA", "20_21delinsCC", None),
    "replacement_unequal": ("ATGAAACGTTCCCCTGCAAGGTTAA", "11_12delinsCCC", None),
    "duplication_CG": ("ATGAAACGCGTTGGCTGCAAGGTTAA", "7_8dup", None),
}


def apply(reference, description):
    """Apply one Synapse HGVS change (c.<...>) to the reference; for the equivalence check."""
    body = description[2:]
    location, _, rest = body.partition("del") if "del" in body else (body, "", "")
    if body.endswith("dup"):
        span = body[:-3].split("_")
        first, last = int(span[0]), int(span[-1])
        return reference[:last] + reference[first - 1:last] + reference[last:]
    if ">" in body:
        position = int(body[:-3])
        return reference[:position - 1] + body[-1] + reference[position:]
    if "ins" in body and "del" not in body:
        left, right = body.split("ins")[0].split("_")
        return reference[:int(left)] + body.split("ins")[1] + reference[int(left):]
    span = location.split("_")
    first, last = int(span[0]), int(span[-1])
    inserted = rest[3:] if rest.startswith("ins") else ""
    return reference[:first - 1] + inserted + reference[last:]


class MutalyzerComparisonTests(unittest.TestCase):
    def test_allele_descriptions_match_mutalyzer(self):
        for name, (observed, mutalyzer, equivalent) in MUTALYZER.items():
            with self.subTest(case=name):
                result = M.analyze_mutations(REFERENCE, observed)
                self.assertEqual(M.hgvs_c_allele(result.mutations, REFERENCE),
                                 "c." + (equivalent or mutalyzer))

    def test_every_description_reproduces_the_observed_sequence(self):
        for name, (observed, _, _) in MUTALYZER.items():
            with self.subTest(case=name):
                description = M.hgvs_c_allele(M.analyze_mutations(REFERENCE, observed).mutations, REFERENCE)
                self.assertEqual(apply(REFERENCE, description), observed)

    def test_insertion_without_reference_uses_both_flanking_positions(self):
        insertion = M.analyze_mutations(REFERENCE, "ATGAAACGTTGGCTGCATTAGGTTAA").mutations[0]
        self.assertEqual(M.mutation_hgvs_c(insertion, REFERENCE), "c.17_18insTT")
        self.assertRegex(M.mutation_hgvs_c(insertion), r"^c\.\d+_\d+insTT$")

    def test_several_changes_form_one_allele(self):
        observed = "ATGAAATGTTGGCTGCAAGCCTAA"
        result = M.analyze_mutations(REFERENCE, observed)
        self.assertEqual(M.hgvs_c_allele(result.mutations, REFERENCE), "c.[7C>T;20_21delinsCC]")
        self.assertEqual(apply(apply(REFERENCE, "c.20_21delinsCC"), "c.7C>T"), observed)

    def test_report_shows_hgvs_descriptions(self):
        result = M.analyze_mutations(REFERENCE, "ATGAAACGTTTTGGCTGCAAGGTTAA")
        report = M.build_mutation_report(result)
        self.assertIn("c.9_10dup", report)
        self.assertNotIn("insTT", report)
        self.assertRegex(report, r"HGVS allele:\s+c\.9_10dup")


if __name__ == "__main__":
    unittest.main()
