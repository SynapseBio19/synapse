"""Known-answer regression tests for the nine P0 findings of the 2026-10-29 audit.

Each test states the scientifically expected behavior. In Stage 0 these tests
were marked ``expectedFailure`` (they reproduced the defects); since Stage 1
they must pass.
"""

from __future__ import annotations

import unittest

from tests._support import PROJECT_ROOT, read_fasta_records, reverse_complement, deterministic_dna

EXAMPLES = dict(read_fasta_records(PROJECT_ROOT / "data" / "examples" / "dna.fasta"))


def _blast_xml(query_from, query_to, align_len, evalue, qseq, hseq, extra_hsp=""):
    return f"""<?xml version="1.0"?>
<BlastOutput><BlastOutput_program>blastn</BlastOutput_program>
<BlastOutput_iterations><Iteration><Iteration_hits><Hit>
<Hit_id>ref|X1|</Hit_id><Hit_def>Test hit</Hit_def><Hit_accession>X1</Hit_accession><Hit_len>500</Hit_len>
<Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>40</Hsp_bit-score><Hsp_evalue>{evalue}</Hsp_evalue>
<Hsp_query-from>{query_from}</Hsp_query-from><Hsp_query-to>{query_to}</Hsp_query-to>
<Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>{align_len}</Hsp_hit-to>
<Hsp_identity>{align_len}</Hsp_identity><Hsp_positive>{align_len}</Hsp_positive><Hsp_gaps>0</Hsp_gaps>
<Hsp_align-len>{align_len}</Hsp_align-len><Hsp_qseq>{qseq}</Hsp_qseq><Hsp_hseq>{hseq}</Hsp_hseq>
<Hsp_midline>{'|' * align_len}</Hsp_midline></Hsp>{extra_hsp}</Hit_hsps></Hit></Iteration_hits></Iteration></BlastOutput_iterations>
</BlastOutput>"""


class P01BlastCoverageAndEvidence(unittest.TestCase):
    def _run(self, xml, query_length):
        from analyses import blast as B
        query = B.build_query("A" * query_length, name="q")
        return B.parse_blast_xml(xml, query=query, rid="T", parameters=B.BlastParameters(),
                                 submitted_at="2026-01-01T00:00:00+00:00")

    def test_query_coverage_is_fraction_of_query(self):
        run = self._run(_blast_xml(1, 32, 32, 1e-5, "A" * 32, "A" * 32), 99)
        self.assertAlmostEqual(run.hits[0].best_coverage_percent, 100 * 32 / 99, places=3)

    def test_non_significant_hit_gets_no_evidence_band(self):
        from analyses import blast as B
        run = self._run(_blast_xml(1, 99, 99, 7.45, "A" * 99, "A" * 99), 99)
        evidence = B.calculate_evidence_score(run.hits[0])
        self.assertEqual(evidence.evidence_band, "NOT_SIGNIFICANT")

    def test_stored_test8_result_recomputed(self):
        # The audited result file: 32-bp HSP on a 99-bp query, E = 7.45.
        import json
        from analyses import blast as B
        data = json.loads((PROJECT_ROOT / "tests" / "fixtures" / "legacy_blast" / "TEST_8_REFERENCE_blast.json").read_text())
        hsp = data["hits"][0]["hsps"][0]
        self.assertEqual(hsp["alignment_length"], 32)
        hit = B.BlastHit(hit_id="h", accession="h", definition="", length=None,
                         hsps=[B.BlastHSP(**{k: hsp[k] for k in B.BlastHSP.__dataclass_fields__ if k in hsp})],
                         query_length=99)
        self.assertAlmostEqual(hit.best_coverage_percent, 32.323, places=2)


class P02CrisprOrfCoordinates(unittest.TestCase):
    def _check(self, strand):
        from analyses import orf as O, crispr as C
        flank = deterministic_dna(40, seed=5)
        cds = "ATG" + "GCTGGACTGAAGTGGCGGGAGAGG" * 3 + "TAA"
        dna = flank + (cds if strand == "+" else reverse_complement(cds)) + deterministic_dna(40, seed=6)
        orf = [o for o in O.find_orfs(dna, include_partial=False, min_length=60)
               if o["strand"] == strand and o["sequence"] == cds][0]
        result = C.run_crispr_analysis(dna, selected_orf=orf)
        self.assertGreater(result.candidate_count, 0)
        for g in result.candidates:
            genomic = dna[g.guide_start - 1:g.guide_end]
            observed = genomic if g.strand == "+" else reverse_complement(genomic)
            self.assertEqual(observed, g.guide, f"{g.candidate_id} on {g.strand}")

    def test_plus_strand_orf(self):
        self._check("+")

    def test_minus_strand_orf(self):
        self._check("-")


def _site(enzyme, flank="AAAA"):
    from analyses import cloning as C
    seq = flank + C.get_recognition_sequence(enzyme) + flank
    analysis = C.analyze_enzyme(seq, enzyme)
    return C.make_site(analysis, analysis["sites"][0])


class P03CloningEndPolarity(unittest.TestCase):
    def test_same_sequence_opposite_polarity_is_incompatible(self):
        from analyses import cloning as C
        self.assertIs(C.ends_compatible(_site("SacI"), _site("HindIII", "CCCC")), False)
        self.assertIs(C.ends_compatible(_site("SphI"), _site("NcoI", "CCCC")), False)

    def test_known_compatible_pairs_remain_compatible(self):
        from analyses import cloning as C
        for a, b in (("XbaI", "SpeI"), ("SalI", "XhoI"), ("BamHI", "BamHI"), ("EcoRV", "SmaI"),
                     ("PstI", "PstI"), ("MspI", "TaqI")):
            self.assertIs(C.ends_compatible(_site(a), _site(b, "CCCC")), True, f"{a}+{b}")


class P04EnzymeData(unittest.TestCase):
    def test_clai_cuts_at_ct_cgat_with_cg_overhang(self):
        from analyses import restriction as R
        site = R.analyze_restriction_site("AAAAATCGATAAAAA", "ClaI", 4)
        self.assertEqual(site["overhang_sequence"], "CG")
        self.assertEqual(site["overhang_length"], 2)

    def test_cloning_overhangs_derive_from_cut_positions(self):
        from analyses import cloning as C
        self.assertEqual(C.get_overhang("NdeI"), "TA")
        self.assertEqual(C.get_overhang("ClaI"), "CG")
        self.assertIs(C.ends_compatible(_site("ClaI"), _site("TaqI", "CCCC")), True)
        self.assertIs(C.ends_compatible(_site("ClaI"), _site("XhoI", "CCCC")), False)


class P05AmbiguousSitesNotCut(unittest.TestCase):
    def test_n_run_produces_no_physical_cut(self):
        from analyses import restriction as R
        digest = R.perform_digest("A" * 50 + "N" * 10 + "A" * 50, ["EcoRI", "BamHI", "HindIII"], "linear")
        self.assertEqual(digest["cut_count"], 0)
        self.assertEqual(len(digest["fragments"]), 1)
        self.assertGreater(len(digest["possible_cuts"]), 0)


class P06StartLoss(unittest.TestCase):
    def test_internal_met_change_is_missense(self):
        from analyses import mutation as M
        ref = "ATG" + "GCT" * 8 + "ATG" + "GCT" * 5 + "TAA"
        alt = "ATG" + "GCT" * 8 + "ATA" + "GCT" * 5 + "TAA"
        c = M.analyze_coding_consequences(ref, alt)
        self.assertEqual((c.missense_count, c.start_loss_count), (1, 0))

    def test_first_codon_change_is_start_loss(self):
        from analyses import mutation as M
        c = M.analyze_coding_consequences("ATGGCTGCTTAA", "ATAGCTGCTTAA")
        self.assertEqual((c.missense_count, c.start_loss_count), (0, 1))


class P07VariantNormalization(unittest.TestCase):
    REF = "GGCAAAAGTC"  # positions 1..10; A-run at 4..7

    def _norm(self, pos, ref, alt):
        from analyses import variant_analysis as V
        v = V.left_normalize_variant(V.Variant("chr1", pos, ref, alt), reference_sequence=self.REF)
        return v.position, v.reference, v.alternative

    def test_non_repeat_deletion_keeps_identity(self):
        self.assertEqual(self._norm(7, "AG", "A"), (7, "AG", "A"))

    def test_insertion_in_run_is_left_shifted_with_anchor(self):
        self.assertEqual(self._norm(7, "A", "AA"), (3, "C", "CA"))

    def test_deletion_in_run_is_left_shifted_with_anchor(self):
        self.assertEqual(self._norm(6, "AA", "A"), (3, "CA", "C"))

    def test_mnv_and_snv_unchanged(self):
        self.assertEqual(self._norm(2, "G", "T"), (2, "G", "T"))
        self.assertEqual(self._norm(4, "AA", "TT"), (4, "AA", "TT"))

    def test_deletion_reaching_sequence_start_is_right_anchored(self):
        # Reference AAAC...: deleting one A of the leading run.
        from analyses import variant_analysis as V
        v = V.left_normalize_variant(V.Variant("chr1", 2, "AA", "A"), reference_sequence="AAACGT")
        self.assertEqual((v.position, v.reference, v.alternative), (1, "AA", "A"))

    def test_reference_mismatch_is_reported(self):
        from analyses import variant_analysis as V
        v = V.left_normalize_variant(V.Variant("chr1", 2, "T", "C"), reference_sequence=self.REF)
        self.assertIn("REF_MISMATCH", v.annotation.get("validation_flags", []))


class P08AlternativeStartTranslation(unittest.TestCase):
    def test_orf_translation_uses_initiator_met(self):
        from analyses import orf as O
        self.assertEqual(O.translate_orf_sequence("GTGAAATTTTAA"), "MKF")
        self.assertEqual(O.translate_orf_sequence("TTGAAATTTTAA"), "MKF")

    def test_protein_translation_from_orf_start(self):
        from analyses import protein as P
        self.assertEqual(P.analyze_translation("GUGAAAUUUUAA", initiator_start=True)["protein"], "MKF")
        # Without the flag (arbitrary sequence, not an ORF start) GUG stays Val.
        self.assertEqual(P.analyze_translation("GUGAAAUUUUAA")["protein"], "VKF")


class P09OrfDefaults(unittest.TestCase):
    def test_reference_cds_ranks_first(self):
        from analyses import orf as O
        orfs = O.find_orfs(EXAMPLES["TEST_8_REFERENCE"])
        self.assertEqual((orfs[0]["start"], orfs[0]["length"], orfs[0]["start_codon"]), (0, 99, "ATG"))

    def test_default_minimum_length_is_75_nt(self):
        from analyses import orf as O
        self.assertEqual(O.DEFAULT_MIN_ORF_LENGTH, 75)
        orfs = O.find_orfs(deterministic_dna(2000, seed=9))
        self.assertTrue(all(o["length"] >= 75 for o in orfs))
        self.assertTrue(all(o["start_codon"] == "ATG" for o in orfs))


# All nine P0 defects were fixed in Stage 1 (2026-10-29); no test is marked
# expectedFailure any more. See docs/CHANGES.md.


if __name__ == "__main__":
    unittest.main()
