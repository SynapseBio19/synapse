"""Phase 6 Tier 1 tests: silent scientific errors (docs/PHASE6_INITIAL_AUDIT.md)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tests import _support  # noqa: F401


# ----------------------------------------------------------------------
# D12: contig names in REF checking and normalization
# ----------------------------------------------------------------------
class VariantContigMatchingTests(unittest.TestCase):
    # A-run at 4..8: deleting one A left-aligns to 3 GA>G (Tan et al. 2015).
    REF = "GGGAAAAAGGG"

    def _deletion(self, chromosome):
        from analyses import variant_analysis as V
        return V.Variant(chromosome, 7, "AA", "A")

    def test_unmatched_contig_is_flagged_not_silently_trimmed(self):
        from analyses import variant_analysis as V
        original = self._deletion("1")
        [result] = V.normalize_variants([original], contig_references={"chr1": self.REF})
        self.assertEqual((result.position, result.reference, result.alternative), (7, "AA", "A"))
        self.assertIn("REFERENCE_CONTIG_NOT_FOUND", result.annotation["validation_flags"])
        self.assertTrue(any("not left-aligned" in note for note in result.normalization_notes))
        self.assertTrue(any("chr1" in note for note in result.normalization_notes))  # suggestion
        # The input variant is not modified by normalization.
        self.assertNotIn("validation_flags", original.annotation)

    def test_reference_check_reports_missing_contig_and_suggestion(self):
        from analyses import variant_analysis as V
        report = V.check_reference_alleles([self._deletion("1")], {"chr1": self.REF})
        self.assertEqual((report["checked"], report["mismatched"], report["not_checked"]), (0, 0, 1))
        self.assertEqual(report["contigs_not_found"], ["1"])
        self.assertEqual(report["suggested_aliases"], {"1": "chr1"})
        self.assertEqual(report["aliased_contigs"], {})

    def test_chr_prefix_mode_matches_both_directions_and_keeps_identity(self):
        from analyses import variant_analysis as V
        for chromosome, reference_name in (("1", "chr1"), ("chr1", "1")):
            with self.subTest(chromosome=chromosome):
                [result] = V.normalize_variants(
                    [self._deletion(chromosome)],
                    contig_references={reference_name: self.REF},
                    contig_matching="chr_prefix",
                )
                self.assertEqual((result.position, result.reference, result.alternative), (3, "GA", "G"))
                self.assertEqual(result.chromosome, chromosome)  # CHROM is never renamed
                self.assertEqual(result.annotation["reference_contig"], reference_name)
                self.assertIn("CONTIG_NAME_ALIASED", result.annotation["validation_flags"])
                self.assertNotIn("REFERENCE_CONTIG_NOT_FOUND", result.annotation["validation_flags"])

    def test_reference_check_counts_aliased_contig(self):
        from analyses import variant_analysis as V
        report = V.check_reference_alleles([self._deletion("1")], {"chr1": self.REF},
                                           contig_matching="chr_prefix")
        self.assertEqual((report["checked"], report["not_checked"]), (1, 0))
        self.assertEqual(report["aliased_contigs"], {"1": "chr1"})

    def test_mitochondrial_names_are_never_aliased_or_suggested(self):
        # hg19 chrM (16,571 bp) and rCRS MT (16,569 bp) are different sequences.
        from analyses import variant_analysis as V
        for chromosome, reference_name in (("MT", "chrMT"), ("chrMT", "MT"), ("M", "chrM"), ("chrM", "M")):
            with self.subTest(chromosome=chromosome):
                report = V.check_reference_alleles(
                    [V.Variant(chromosome, 7, "AA", "A")], {reference_name: self.REF},
                    contig_matching="chr_prefix",
                )
                self.assertEqual(report["contigs_not_found"], [chromosome])
                self.assertEqual(report["aliased_contigs"], {})
                self.assertEqual(report["suggested_aliases"], {})

    def test_case_variants_are_not_folded(self):
        from analyses import variant_analysis as V
        report = V.check_reference_alleles([self._deletion("Chr1")], {"chr1": self.REF, "1": self.REF},
                                           contig_matching="chr_prefix")
        self.assertEqual(report["contigs_not_found"], ["Chr1"])

    def test_exact_name_is_preferred_over_alias(self):
        from analyses import variant_analysis as V
        references = {"1": self.REF, "chr1": "CCCCCCCCCCC"}
        [result] = V.normalize_variants([self._deletion("1")], contig_references=references,
                                        contig_matching="chr_prefix")
        self.assertEqual(result.position, 3)
        self.assertNotIn("CONTIG_NAME_ALIASED", result.annotation.get("validation_flags", []))

    def test_unknown_matching_mode_is_rejected(self):
        from analyses import variant_analysis as V
        with self.assertRaises(V.VariantAnalysisError):
            V.normalize_variants([self._deletion("1")], contig_references={"1": self.REF},
                                 contig_matching="fuzzy")

    def test_api_reports_unnormalized_variants_truthfully(self):
        from core import api
        vcf = ("##fileformat=VCFv4.3\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
               "1\t7\t.\tAA\tA\t.\tPASS\t.\n")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "v.vcf"
            path.write_text(vcf, encoding="utf-8")
            exact = api.analyze_vcf(str(path), contig_references={"chr1": self.REF})
            aliased = api.analyze_vcf(str(path), contig_references={"chr1": self.REF},
                                      contig_matching="chr_prefix")

        self.assertFalse(exact.parameters["normalized"])
        self.assertEqual(exact.parameters["contig_matching"], "exact")
        self.assertEqual(exact.data["normalization"]["variants_without_reference"], 1)
        self.assertEqual(exact.data["normalization"]["contigs_not_found"], ["1"])
        self.assertTrue(any("'1'" in w and "chr1" in w for w in exact.warnings))
        self.assertEqual(exact.data["variants"][0].position, 7)

        self.assertTrue(aliased.parameters["normalized"])
        self.assertEqual(aliased.data["normalization"]["variants_without_reference"], 0)
        self.assertEqual(aliased.data["normalization"]["aliased_contigs"], {"1": "chr1"})
        self.assertTrue(any("prefix" in w for w in aliased.warnings))
        self.assertEqual(aliased.data["variants"][0].position, 3)


if __name__ == "__main__":
    unittest.main()
