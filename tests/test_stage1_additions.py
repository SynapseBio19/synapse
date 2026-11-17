"""Additional tests for Stage 1 behavior (ORF policy, ORF adapter, BLAST coverage)."""

from __future__ import annotations

import json
import unittest

from tests import _support
from tests._support import deterministic_dna, reverse_complement
from tests.golden_scenarios import EXAMPLES, RANDOM_600, _orf_view


class OrfLegacyModeTests(unittest.TestCase):
    """The pre-2026-10-29 ORF behavior must stay reproducible on request."""

    def test_legacy_settings_reproduce_baseline_snapshot(self):
        from analyses import orf as O
        baseline = _support.load_golden("orf_legacy_baseline")
        for name in ("TEST_8_REFERENCE", "TEST_8_REVERSE", "TEST_8_SNP_NONSENSE"):
            observed = json.loads(json.dumps(_support.normalize_for_json(
                [_orf_view(o) for o in O.find_orfs(EXAMPLES[name], **O.LEGACY_ORF_SETTINGS)])))
            self.assertEqual(baseline[name], observed, name)
        result = O.analyze_orfs_enhanced(RANDOM_600, **O.LEGACY_ORF_SETTINGS)
        self.assertEqual(baseline["random600_count"], result["candidate_count"])

    def test_alternative_start_policy_finds_gtg_orf(self):
        from analyses import orf as O
        seq = "CC" + "GTG" + "GCT" * 30 + "TAA" + "CC"
        self.assertEqual(O.find_orfs(seq, include_partial=False), [])
        orfs = O.find_orfs(seq, include_partial=False, start_codon_policy="atg_and_alternative")
        self.assertEqual(orfs[0]["start_codon"], "GTG")
        self.assertTrue(O.translate_orf_sequence(orfs[0]["sequence"]).startswith("MA"))

    def test_nested_flag(self):
        from analyses import orf as O
        seq = "ATG" + "GCT" * 10 + "ATG" + "GCT" * 30 + "TAA"
        self.assertEqual(len(O.find_orfs(seq, include_partial=False)), 1)
        self.assertEqual(len(O.find_orfs(seq, include_partial=False, nested=True)), 2)

    def test_parameters_recorded(self):
        from analyses import orf as O
        result = O.analyze_orfs(deterministic_dna(300))
        params = result["parameters"]
        self.assertEqual(params["min_length_nt"], 75)
        self.assertEqual(params["start_codons"], ["ATG"])
        self.assertFalse(params["nested"])

    def test_unknown_policy_rejected(self):
        from analyses import orf as O
        with self.assertRaises(ValueError):
            O.find_orfs("ATGAAATAA", start_codon_policy="anything")


class OrfAdapterTests(unittest.TestCase):
    def test_synapse_orf_is_zero_based(self):
        from analyses import orf as O
        from core.orf_context import resolve_orf_region
        dna = "GG" + "ATG" + "GCT" * 30 + "TAA" + "GG"
        orf = O.find_orfs(dna, include_partial=False)[0]
        region = resolve_orf_region(orf, dna)
        assert region is not None
        self.assertTrue(region.validated)
        self.assertEqual((region.start_1based, region.end_1based, region.strand), (3, 98, "+"))

    def test_minus_strand_region_is_forward_genomic(self):
        from analyses import orf as O
        from core.orf_context import resolve_orf_region
        cds = "ATG" + "GCT" * 30 + "TAA"
        dna = "GG" + reverse_complement(cds) + "GG"
        orf = [o for o in O.find_orfs(dna, include_partial=False) if o["strand"] == "-"][0]
        region = resolve_orf_region(orf, dna)
        assert region is not None
        self.assertTrue(region.validated)
        self.assertEqual(region.forward_region, dna[2:98])
        self.assertEqual(region.coding_sequence, cds)

    def test_generic_dict_keeps_one_based_start(self):
        from core.orf_context import resolve_orf_region
        region = resolve_orf_region({"sequence": "ATGCGT", "start": 101})
        assert region is not None
        self.assertEqual((region.start_1based, region.end_1based), (101, 106))

    def test_mismatching_parent_is_not_trusted(self):
        from core.orf_context import resolve_orf_region
        orf = {"sequence": "ATGAAA", "start": 0, "stop": 5, "strand": "+", "frame": 1}
        region = resolve_orf_region(orf, "CCCCCCCC")
        assert region is not None and region.note is not None
        self.assertFalse(region.validated)
        self.assertIn("does not match", region.note)

    def test_blast_orf_query_region(self):
        from analyses import orf as O, blast as B
        cds = "ATG" + "GCT" * 30 + "TAA"
        dna = "GGGGG" + reverse_complement(cds) + "GG"
        orf = [o for o in O.find_orfs(dna, include_partial=False) if o["strand"] == "-"][0]
        query = B.build_orf_query(orf, parent_sequence=dna)
        self.assertEqual(query.sequence, cds)
        self.assertEqual((query.region_start, query.region_end), (6, 101))
        assert query.feature_name is not None
        self.assertIn("- strand", query.feature_name)


class BlastCoverageTests(unittest.TestCase):
    def _hsp(self, qfrom, qto, evalue=1e-30):
        from analyses import blast as B
        n = qto - qfrom + 1
        return B.BlastHSP(number=1, bit_score=80.0, raw_score=None, evalue=evalue,
                          query_start=qfrom, query_end=qto, subject_start=1, subject_end=n,
                          query_frame=1, subject_frame=1, identity=n, positive=n, gaps=0,
                          alignment_length=n, query_sequence="A" * n, subject_sequence="A" * n,
                          midline="|" * n)

    def test_union_of_hsps(self):
        from analyses import blast as B
        hit = B.BlastHit("h", "h", "", 1000, hsps=[self._hsp(1, 40), self._hsp(31, 60), self._hsp(81, 100)],
                         query_length=100)
        self.assertAlmostEqual(hit.best_coverage_percent, 80.0)
        self.assertAlmostEqual(hit.hsps[0].query_coverage, 40.0)

    def test_minus_strand_hsp_interval(self):
        hsp = self._hsp(1, 50)
        hsp.query_start, hsp.query_end = 50, 1
        hsp.query_length = 100
        self.assertAlmostEqual(hsp.query_coverage, 50.0)

    def test_threshold_is_configurable(self):
        from analyses import blast as B
        hit = B.BlastHit("h", "h", "", 1000, hsps=[self._hsp(1, 100, evalue=0.01)], query_length=100)
        self.assertEqual(B.calculate_evidence_score(hit).evidence_band, "NOT_SIGNIFICANT")
        relaxed = B.calculate_evidence_score(hit, significance_evalue=0.05)
        self.assertNotEqual(relaxed.evidence_band, "NOT_SIGNIFICANT")


class RestrictionPolarityTests(unittest.TestCase):
    def test_restriction_module_polarity(self):
        from analyses import restriction as R
        sac = R.analyze_restriction_site("AAGAGCTCAA", "SacI", 2)
        hind = R.analyze_restriction_site("AAAAGCTTAA", "HindIII", 2)
        self.assertFalse(R.sticky_ends_are_compatible(sac, hind))

    def test_possible_cuts_can_be_included(self):
        from analyses import restriction as R
        d = R.perform_digest("A" * 50 + "N" * 10 + "A" * 50, ["EcoRI"], "linear", include_possible=True)
        self.assertGreater(d["cut_count"], 0)
        self.assertTrue(all(c["status"] == "POSSIBLE" for c in d["cuts"]))


if __name__ == "__main__":
    unittest.main()
