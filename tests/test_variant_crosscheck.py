"""P2-4: designed primers and guides checked against variants.

A variant under a primer reduces binding, most of all at the 3' end where
the polymerase extends; a variant in a guide's seed or PAM most reduces
cleavage. The check reports every overlap and marks those critical parts.
"""

from __future__ import annotations

import unittest

from tests import _support  # noqa: F401
from tests._support import deterministic_dna

DNA = deterministic_dna(900)


def _variant(position, reference, alternative="A", chromosome="seq"):
    from analyses.variant_analysis import Variant
    return Variant(chromosome, position, reference, alternative)


class PrimerCrossCheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from core import api
        cls.result = api.design_primers(DNA)
        pair = cls.result.data.ranked_pairs[0].pair
        cls.forward, cls.reverse = pair.forward, pair.reverse

    def _check(self, *variants):
        from core import api
        return api.cross_check_variants(self.result, list(variants), top_n=1)

    def test_binding_sites_follow_primer_orientation(self):
        from analyses.primer import primer_pair_binding_sites
        forward, reverse = primer_pair_binding_sites(self.result.data.ranked_pairs[0])
        self.assertEqual((forward["start"], forward["end"]), (self.forward.start + 1, self.forward.end))
        self.assertEqual(forward["critical"], (self.forward.end - 4, self.forward.end))
        self.assertEqual(reverse["critical"], (self.reverse.start + 1, self.reverse.start + 5))

    def test_a_variant_at_the_3_prime_end_is_high(self):
        last = self.forward.end
        (finding,) = self._check(_variant(last, DNA[last - 1]))
        self.assertEqual((finding["feature"], finding["region"], finding["severity"]),
                         ("pair 1 forward primer", "3' end", "HIGH"))
        first_of_reverse = self.reverse.start + 1
        (finding,) = self._check(_variant(first_of_reverse, DNA[first_of_reverse - 1]))
        self.assertEqual((finding["feature"], finding["severity"]), ("pair 1 reverse primer", "HIGH"))

    def test_a_variant_in_the_5_prime_part_is_moderate(self):
        position = self.forward.start + 1
        (finding,) = self._check(_variant(position, DNA[position - 1]))
        self.assertEqual((finding["region"], finding["severity"]), ("binding site", "MODERATE"))

    def test_a_deletion_reaching_into_the_site_counts(self):
        position = self.forward.start - 2                     # REF covers 3 bases into the primer
        reference = DNA[position - 1:position + 5]
        (finding,) = self._check(_variant(position, reference))
        self.assertEqual(finding["feature"], "pair 1 forward primer")

    def test_outside_and_other_contigs_are_not_reported(self):
        outside = self.forward.start - 30
        self.assertEqual(self._check(_variant(outside, DNA[outside - 1])), [])
        from core import api
        inside = _variant(self.forward.end, DNA[self.forward.end - 1], chromosome="other")
        self.assertEqual(api.cross_check_variants(self.result, [inside], chromosome="seq", top_n=1), [])


class GuideCrossCheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from core import api
        cls.result = api.design_crispr_guides(DNA)
        cls.plus = next(c for c in cls.result.data.candidates if c.strand == "+")
        cls.minus = next(c for c in cls.result.data.candidates if c.strand == "-")

    def _findings(self, candidate, position):
        from analyses.crispr import guide_binding_site
        from analyses.variant_analysis import check_features_against_variants
        return check_features_against_variants([guide_binding_site(candidate)],
                                               [_variant(position, DNA[position - 1])])

    def test_pam_seed_and_distal_part_on_the_plus_strand(self):
        guide = self.plus                      # PAM to the right of the protospacer
        self.assertEqual(self._findings(guide, guide.pam_end)[0]["severity"], "HIGH")
        self.assertEqual(self._findings(guide, guide.guide_end)[0]["region"], "seed or PAM")
        self.assertEqual(self._findings(guide, guide.guide_start)[0]["severity"], "MODERATE")

    def test_pam_seed_and_distal_part_on_the_minus_strand(self):
        guide = self.minus                     # PAM to the left in sequence coordinates
        self.assertLess(guide.pam_end, guide.guide_start)
        self.assertEqual(self._findings(guide, guide.pam_start)[0]["severity"], "HIGH")
        self.assertEqual(self._findings(guide, guide.guide_start)[0]["severity"], "HIGH")
        self.assertEqual(self._findings(guide, guide.guide_end)[0]["severity"], "MODERATE")

    def test_the_api_checks_the_ranked_guides(self):
        from core import api
        guide = self.result.data.candidates[0]
        findings = api.cross_check_variants(self.result, [_variant(guide.guide_start, DNA[guide.guide_start - 1])],
                                            top_n=1)
        self.assertEqual([f["feature"] for f in findings], [f"guide {guide.candidate_id}"])

    def test_other_results_are_refused(self):
        from core import api
        with self.assertRaises(ValueError):
            api.cross_check_variants(api.find_orfs(DNA), [])


if __name__ == "__main__":
    unittest.main()
