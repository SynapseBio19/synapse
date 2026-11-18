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


class MenuIntegrationTests(unittest.TestCase):
    """T12: the dialogs check designs against the variants mapped onto the sequence."""

    def _drive(self, answers, call):
        import builtins
        import contextlib
        import io
        from unittest import mock
        pending = list(answers)
        buffer = io.StringIO()
        with mock.patch.object(builtins, "input", lambda prompt="": pending.pop(0)), \
                contextlib.redirect_stdout(buffer):
            result = call()
        return result, buffer.getvalue()

    def test_primer_dialog_reports_a_variant_at_a_3_prime_end(self):
        from core import api
        from core.pipeline import SynapsePipeline
        from ui import primer_menu
        forward = api.design_primers(DNA).data.ranked_pairs[0].pair.forward
        pipeline = SynapsePipeline(DNA, "seq")
        last_base = forward.end            # 1-based position of the forward primer's 3' base
        pipeline.context.set_mapped_variants("seq", [
            _variant(last_base, DNA[last_base - 1], "A" if DNA[last_base - 1] != "A" else "C"),
            _variant(5, DNA[4], "G", chromosome="other"),   # another contig: ignored
        ])
        _, output = self._drive(["", "", "", "", ""],      # Tm model, four conditions: defaults
                                lambda: primer_menu._primer_ui_design(DNA, 0, len(DNA), pipeline))
        self.assertIn("VARIANTS IN BINDING SITES", output)
        self.assertIn("Variants on the contig:", output)
        line = next(text for text in output.splitlines() if text.startswith("pair 1 forward primer"))
        self.assertIn(f"{last_base} ", line)
        self.assertIn("HIGH", line)
        self.assertIn("3' end", line)

    def test_no_check_without_variants_or_on_another_sequence(self):
        from core.pipeline import SynapsePipeline
        from ui import primer_menu
        pipeline = SynapsePipeline(DNA, "seq")
        _, output = self._drive(["", "", "", "", ""],
                                lambda: primer_menu._primer_ui_design(DNA, 0, len(DNA), pipeline))
        self.assertNotIn("VARIANTS IN BINDING SITES", output)
        pipeline.context.set_mapped_variants("seq", [_variant(100, DNA[99])])
        other = deterministic_dna(900, seed=5)
        self.assertIsNone(pipeline.context.variants_for_sequence(other))
        self.assertIsNotNone(pipeline.context.variants_for_sequence(DNA.lower()))

    def test_guide_design_is_checked_through_the_shared_section(self):
        import contextlib
        import io
        from analyses.crispr import design_guides
        from core.pipeline import SynapsePipeline
        from ui.variant_check import print_variant_check
        result = design_guides(DNA)
        guide = result.candidates[0]
        pipeline = SynapsePipeline(DNA, "seq")
        pipeline.context.set_mapped_variants("seq", [_variant(guide.pam_start, DNA[guide.pam_start - 1])])
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            findings = print_variant_check(pipeline, DNA, "crispr", result)
        self.assertTrue(findings)
        assert findings is not None
        self.assertEqual(findings[0]["severity"], "HIGH")
        self.assertIn("VARIANTS IN BINDING SITES", buffer.getvalue())

    def test_crispr_menu_accepts_the_pipeline(self):
        import inspect
        from ui import crispr_menu
        self.assertIn("pipeline", inspect.signature(crispr_menu.crispr_menu).parameters)


class VariantDialogMappingTests(unittest.TestCase):
    def test_normalizing_with_a_mapped_contig_keeps_the_variants_for_the_sequence(self):
        import builtins
        import contextlib
        import io
        from unittest import mock
        from core.pipeline import SynapsePipeline
        from ui import variant_analysis_menu
        pipeline = SynapsePipeline(DNA, "seq")
        records = [f"seq\t{p}\t.\t{DNA[p - 1]}\t{'A' if DNA[p - 1] != 'A' else 'C'}\t50\tPASS\t." for p in (50, 400)]
        vcf = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"] + records
        # load: paste text, empty line ends it; normalize, map "seq"; leave
        answers = ["1", ""] + vcf + ["", "3", "seq", "0"]
        with mock.patch.object(builtins, "input", lambda prompt="": answers.pop(0)), \
                contextlib.redirect_stdout(io.StringIO()) as out:
            variant_analysis_menu.run_variant_analysis_ui(pipeline=pipeline)
        mapping = pipeline.context.variants_for_sequence(DNA)
        self.assertIsNotNone(mapping, out.getvalue()[-1500:])
        assert mapping is not None
        self.assertEqual(mapping["contig"], "seq")
        self.assertEqual(sorted(v.position for v in mapping["variants"]), [50, 400])
        self.assertIn("later primer and guide designs on it are checked", out.getvalue())
