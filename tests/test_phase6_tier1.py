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


def _filler(length):
    # A/C only: contains no recognition site of any enzyme in the database.
    return ("CA" * length)[:length]


# ----------------------------------------------------------------------
# D4: cloning uses the shared matcher; ambiguity never makes a site
# ----------------------------------------------------------------------
class CloningAmbiguousSiteTests(unittest.TestCase):
    N_RUN = "A" * 10 + "N" * 6 + "A" * 10

    def test_n_run_is_possible_not_confirmed(self):
        from analyses import cloning as C
        self.assertEqual(C.find_restriction_sites(self.N_RUN, "GAATTC"), [])
        self.assertEqual(C.find_restriction_sites(self.N_RUN, "GAATTC", include_possible=True), [11])
        analysis = C.analyze_enzyme(self.N_RUN, "EcoRI")
        self.assertEqual((analysis["site_count"], analysis["sites"]), (0, []))
        self.assertEqual((analysis["possible_site_count"], analysis["possible_sites"]), (1, [11]))

    def test_cloning_and_restriction_report_identical_sites(self):
        from analyses import cloning as C, restriction as R
        from tests._support import deterministic_dna
        sequence = deterministic_dna(2000, seed=99)
        sequence = sequence[:500] + "N" + sequence[501:1200] + "RY" + sequence[1202:]
        for enzyme in sorted(C.RESTRICTION_ENZYMES):
            matches = R.find_restriction_site_matches(sequence, enzyme)
            analysis = C.analyze_enzyme(sequence, enzyme)
            with self.subTest(enzyme=enzyme):
                self.assertEqual(analysis["sites"],
                                 [m["position_1based"] for m in matches if m["status"] == "CONFIRMED"])
                self.assertEqual(analysis["possible_sites"],
                                 [m["position_1based"] for m in matches if m["status"] == "POSSIBLE"])

    def test_cloning_has_no_private_recognition_matcher(self):
        from analyses import cloning as C, restriction as R
        self.assertFalse(hasattr(C, "_recognition_matches"))
        self.assertIs(C.bases_match, R.bases_match)

    def _vector_with_phantom_ecori(self):
        # The only EcoRI-compatible site in the vector depends on an N.
        return _filler(40) + "GANTTC" + _filler(40) + "AAGCTT" + _filler(40)

    INSERT = _filler(10) + "GAATTC" + _filler(60) + "AAGCTT" + _filler(10)

    def test_ambiguity_dependent_site_does_not_define_a_vector_end(self):
        from analyses import cloning as C
        result = C.evaluate_cloning_strategy(self._vector_with_phantom_ecori(), self.INSERT,
                                             "EcoRI", "HindIII", "EcoRI", "HindIII")
        self.assertEqual(result["status"], "INVALID_VECTOR_ENDS")
        self.assertIn("ambiguity", result["reason"])

    def test_possible_site_prevents_a_uniqueness_claim(self):
        # One confirmed EcoRI site plus one ambiguity-dependent site: the
        # enzyme may cut twice, so it cannot define a single clean end.
        from analyses import cloning as C
        vector = (_filler(30) + "GAATTC" + _filler(30) + "GANTTC" + _filler(30)
                  + "AAGCTT" + _filler(30))
        result = C.evaluate_cloning_strategy(vector, self.INSERT, "EcoRI", "HindIII", "EcoRI", "HindIII")
        self.assertEqual(result["status"], "INVALID_VECTOR_ENDS")
        self.assertIn("ambiguity", result["reason"])

    def test_single_enzyme_and_suitability_ignore_phantom_sites(self):
        from analyses import cloning as C
        vector = _filler(40) + "GANTTC" + _filler(40)
        insert = _filler(10) + "GAATTC" + _filler(10)
        single = C.analyze_cloning(vector, insert, "EcoRI", "EcoRI")
        self.assertFalse(single["validity"]["valid"])
        self.assertEqual(single["vector"]["site_count"], 0)
        self.assertTrue(any("ambiguity" in w for w in single["warnings"]))
        suitability = C.assess_enzyme_suitability(vector, insert, "EcoRI")
        self.assertFalse(suitability["suitable"])
        self.assertTrue(any("ambiguity" in r for r in suitability["reasons"]))

    def test_planner_does_not_build_a_strategy_on_a_phantom_site(self):
        from analyses import cloning as C
        planner = C.automatic_cloning_planner(self._vector_with_phantom_ecori(), self.INSERT,
                                              enzyme_names=["EcoRI", "HindIII"], top_n=5)
        self.assertEqual(planner["valid_count"], 0)
        diagnostics = C._site_diagnostics(self._vector_with_phantom_ecori(), ["EcoRI"], "vector")
        self.assertEqual(diagnostics[0]["status"], "AMBIGUITY_DEPENDENT")
        self.assertEqual(diagnostics[0]["possible_site_count"], 1)


if __name__ == "__main__":
    unittest.main()
