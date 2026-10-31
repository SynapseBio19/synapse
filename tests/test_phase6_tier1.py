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


# ----------------------------------------------------------------------
# D3: restriction map and site-based results honor circular topology
# ----------------------------------------------------------------------
class RestrictionCircularTopologyTests(unittest.TestCase):
    # EcoRI G^AATTC split across the origin: ...GAA | TTC...
    ECORI_ACROSS = "TTC" + _filler(120) + "GAA"
    # BamHI G^GATCC whose top-strand cut falls exactly at the origin.
    BAMHI_AT_ORIGIN = "GATCC" + _filler(120) + "G"

    def test_map_reports_origin_spanning_site_only_when_circular(self):
        from analyses import restriction as R
        length = len(self.ECORI_ACROSS)
        self.assertEqual(R.create_restriction_map(self.ECORI_ACROSS, ["EcoRI"]), [])  # linear: no site
        [site] = R.create_restriction_map(self.ECORI_ACROSS, ["EcoRI"], topology="circular")
        self.assertTrue(site["wraps_origin"])
        self.assertEqual(site["match_status"], "CONFIRMED")
        self.assertEqual(site["site_position"], length - 2)          # 1-based start of GAATTC
        self.assertEqual(site["site_end_1based"], 3)
        self.assertEqual((site["top_cut_0based"], site["bottom_cut_0based"]), (length - 2, 2))
        self.assertEqual((site["top_cut"], site["bottom_cut"]), (length - 1, 3))
        self.assertEqual(site["overhang_sequence"], "AATT")

    def test_cut_exactly_at_the_origin(self):
        from analyses import restriction as R
        [site] = R.create_restriction_map(self.BAMHI_AT_ORIGIN, ["BamHI"], topology="circular")
        self.assertEqual((site["top_cut_0based"], site["bottom_cut_0based"]), (0, 4))
        self.assertEqual(site["overhang_sequence"], "GATC")
        self.assertEqual(R.perform_digest(self.BAMHI_AT_ORIGIN, ["BamHI"], "circular")["cut_positions"], [0])

    def test_analysis_agrees_with_its_own_digest(self):
        from analyses import restriction as R
        circular = R.analyze_restriction(self.ECORI_ACROSS, ["EcoRI", "BamHI"], "circular")
        self.assertEqual(circular["digest"]["cut_count"], 1)
        self.assertEqual(len(circular["restriction_map"]), 1)
        self.assertEqual(circular["enzyme_results"]["EcoRI"]["confirmed_site_count"], 1)
        self.assertEqual(circular["unique_cutters"], ["EcoRI"])
        self.assertEqual(circular["non_cutters"], ["BamHI"])
        linear = R.analyze_restriction(self.ECORI_ACROSS, ["EcoRI", "BamHI"], "linear")
        self.assertEqual(linear["digest"]["cut_count"], 0)
        self.assertEqual(linear["non_cutters"], ["EcoRI", "BamHI"])

    def test_enzyme_pair_needs_the_origin_spanning_site(self):
        from analyses import restriction as R
        sequence = "TTC" + _filler(60) + "AAGCTT" + _filler(60) + "GAA"
        pairs = {(p["enzyme_a"], p["enzyme_b"])
                 for p in R.recommend_enzyme_pairs(sequence, ["EcoRI", "HindIII"], topology="circular")}
        self.assertEqual(pairs, {("EcoRI", "HindIII")})
        self.assertEqual(R.recommend_enzyme_pairs(sequence, ["EcoRI", "HindIII"]), [])

    def test_map_and_cutter_lists_match_the_digest_for_every_enzyme(self):
        from analyses import restriction as R
        from tests._support import deterministic_dna
        base = deterministic_dna(1500, seed=21)
        for enzyme in sorted(R.RESTRICTION_ENZYMES):
            sequences = [base]
            confirmed = [m for m in R.find_restriction_site_matches(base, enzyme) if m["status"] == "CONFIRMED"]
            if confirmed:  # rotate so the first site spans the origin
                split = confirmed[0]["position"] + 1
                sequences.append(base[split:] + base[:split])
            for sequence in sequences:
                for topology in ("linear", "circular"):
                    with self.subTest(enzyme=enzyme, topology=topology, rotated=sequence is not base):
                        mapped = sorted({s["top_cut_0based"]
                                         for s in R.create_restriction_map(sequence, [enzyme], topology=topology)
                                         if s["match_status"] == "CONFIRMED"
                                         and (topology == "circular" or 0 < s["top_cut_0based"] < len(sequence))})
                        self.assertEqual(mapped, R.perform_digest(sequence, [enzyme], topology)["cut_positions"])
                        has_site = bool(R.find_restriction_site_matches(sequence, enzyme, topology=topology))
                        self.assertEqual(enzyme in R.find_non_cutters(sequence, [enzyme], topology=topology),
                                         not has_site)

    def test_wrappers_and_cloning_handoff_pass_topology(self):
        from analyses import restriction as R
        self.assertEqual(len(R.get_restriction_map(self.ECORI_ACROSS, ["EcoRI"], topology="circular")), 1)
        self.assertEqual(R.get_unique_cutters(self.ECORI_ACROSS, ["EcoRI"], topology="circular"), ["EcoRI"])
        self.assertEqual(R.get_non_cutters(self.ECORI_ACROSS, ["EcoRI"], topology="circular"), [])
        self.assertEqual(R.summarize_enzyme_sites(self.ECORI_ACROSS, "EcoRI", topology="circular")
                         ["confirmed_positions"], [len(self.ECORI_ACROSS) - 2])
        insert = _filler(10) + "GAATTC" + _filler(10)
        self.assertTrue(R.get_cloning_compatible_sites(self.ECORI_ACROSS, insert, "EcoRI",
                                                       vector_topology="circular")["ready_for_cloning"])
        self.assertFalse(R.get_cloning_compatible_sites(self.ECORI_ACROSS, insert, "EcoRI")["ready_for_cloning"])

    def test_map_menu_uses_the_sequence_topology(self):
        import contextlib
        import io
        from unittest import mock
        from analyses import restriction as R
        answers = iter(["2", "EcoRI", "", "", "0"])  # map, enzyme, keep topology, pause, back
        with mock.patch("builtins.input", lambda *_: next(answers)), \
                mock.patch.object(R, "get_restriction_map", wraps=R.get_restriction_map) as spy, \
                contextlib.redirect_stdout(io.StringIO()):
            R.run_restriction_analysis_ui(self.ECORI_ACROSS, topology="circular")
        self.assertEqual(spy.call_args.kwargs.get("topology"), "circular")


# ----------------------------------------------------------------------
# D5: cloning analyses honor the topology of vector and insert
# ----------------------------------------------------------------------
class CloningCircularTopologyTests(unittest.TestCase):
    # Circular vector: EcoRI split across the origin, HindIII inside.
    VECTOR = "TTC" + _filler(100) + "AAGCTT" + _filler(100) + "GAA"
    # Linear insert ordered like the vector sites: HindIII ... EcoRI.
    INSERT = _filler(10) + "AAGCTT" + _filler(60) + "GAATTC" + _filler(10)
    ENDS = ("HindIII", "EcoRI", "HindIII", "EcoRI")

    def test_site_finding_honors_topology(self):
        from analyses import cloning as C
        length = len(self.VECTOR)
        self.assertEqual(C.find_restriction_sites(self.VECTOR, "GAATTC"), [])
        self.assertEqual(C.find_restriction_sites(self.VECTOR, "GAATTC", topology="circular"), [length - 2])
        analysis = C.analyze_enzyme(self.VECTOR, "EcoRI", topology="circular")
        self.assertEqual((analysis["sites"], analysis["topology"]), ([length - 2], "circular"))
        self.assertEqual(C.analyze_enzyme(self.VECTOR, "EcoRI")["topology"], "linear")

    def test_strategy_uses_the_origin_spanning_vector_site(self):
        from analyses import cloning as C
        linear = C.evaluate_cloning_strategy(self.VECTOR, self.INSERT, *self.ENDS)
        self.assertEqual(linear["status"], "INVALID_VECTOR_ENDS")
        circular = C.evaluate_cloning_strategy(self.VECTOR, self.INSERT, *self.ENDS, vector_topology="circular")
        self.assertEqual(circular["status"], "VALID")
        self.assertEqual((circular["vector"]["topology"], circular["insert"]["topology"]), ("circular", "linear"))
        self.assertTrue(circular["strategy"]["directional"])
        digest = circular["digest"]["vector"]
        self.assertEqual((digest["topology"], digest["fragment_count"]), ("circular", 2))
        self.assertTrue(any("circular" in note for note in circular["assumptions"]))
        self.assertEqual(C.analyze_enzyme_pair(self.VECTOR, self.INSERT, *self.ENDS,
                                               vector_topology="circular")["status"], "VALID")

    def test_circular_insert_donor(self):
        from analyses import cloning as C
        donor = "TTC" + _filler(20) + "AAGCTT" + _filler(20) + "GAA"
        vector = _filler(30) + "AAGCTT" + _filler(30) + "GAATTC" + _filler(30)
        self.assertEqual(C.evaluate_cloning_strategy(vector, donor, *self.ENDS)["status"], "INVALID_INSERT_ENDS")
        self.assertEqual(C.evaluate_cloning_strategy(vector, donor, *self.ENDS,
                                                     insert_topology="circular")["status"], "VALID")

    def test_planner_single_enzyme_and_suitability(self):
        from analyses import cloning as C
        names = ["EcoRI", "HindIII"]
        self.assertEqual(C.automatic_cloning_planner(self.VECTOR, self.INSERT, enzyme_names=names)["valid_count"], 0)
        planner = C.automatic_cloning_planner(self.VECTOR, self.INSERT, enzyme_names=names,
                                              vector_topology="circular")
        self.assertGreater(planner["valid_count"], 0)
        self.assertEqual(planner["topology"], {"vector": "circular", "insert": "linear"})
        self.assertEqual(C.run_cloning_analysis(self.VECTOR, self.INSERT, enzyme_names=names,
                                                vector_topology="circular")["status"], planner["status"])

        vector = "TTC" + _filler(50) + "GAA"
        insert = _filler(10) + "GAATTC" + _filler(10)
        self.assertFalse(C.analyze_cloning(vector, insert, "EcoRI", "EcoRI")["validity"]["valid"])
        self.assertTrue(C.analyze_cloning(vector, insert, "EcoRI", "EcoRI",
                                          vector_topology="circular")["validity"]["valid"])
        self.assertTrue(C.assess_enzyme_suitability(vector, insert, "EcoRI", vector_topology="circular")["suitable"])

    def test_unknown_topology_is_rejected(self):
        from analyses import cloning as C
        with self.assertRaises(ValueError):
            C.evaluate_cloning_strategy(self.VECTOR, self.INSERT, *self.ENDS, vector_topology="supercoiled")

    def test_pipeline_passes_the_loaded_topology(self):
        from core.pipeline import SynapsePipeline
        pipeline = SynapsePipeline(sequence=self.VECTOR, header="vector", topology="circular")
        self.assertEqual(pipeline.run_cloning(self.INSERT, *self.ENDS)["status"], "VALID")
        self.assertGreater(pipeline.run_cloning_planner(self.INSERT, ["EcoRI", "HindIII"])["valid_count"], 0)

    def test_menu_passes_vector_and_insert_topology(self):
        import contextlib
        import io
        from unittest import mock
        from analyses import cloning as C, restriction as R
        real_planner = C.automatic_cloning_planner

        def two_enzyme_planner(vector, insert, **kwargs):
            kwargs["enzyme_names"] = ["EcoRI", "HindIII"]
            return real_planner(vector, insert, **kwargs)

        records = [{"header": "vector", "sequence": self.VECTOR, "topology": "circular"},
                   {"header": "insert", "sequence": self.INSERT, "topology": "linear"}]
        # cloning, insert #2, planner, best strategy, pause, back, back
        answers = iter(["2", "2", "1", "1", "", "0", "0"])
        with mock.patch("builtins.input", lambda *_: next(answers)), \
                mock.patch.object(C, "automatic_cloning_planner", side_effect=two_enzyme_planner) as spy, \
                contextlib.redirect_stdout(io.StringIO()):
            R.run_restriction_cloning_ui(self.VECTOR, "vector", records, topology="circular")
        kwargs = spy.call_args.kwargs
        self.assertEqual((kwargs["vector_topology"], kwargs["insert_topology"]), ("circular", "linear"))


# ----------------------------------------------------------------------
# D1: one initiation-codon rule, taken from the genetic code
# ----------------------------------------------------------------------
class OrfInitiationCodonTests(unittest.TestCase):
    TABLES = (1, 2, 4, 11)

    def test_table1_ctg_initiation_is_methionine_on_every_path(self):
        from analyses import orf as O, protein as P
        from core import api
        sequence = "CTG" + "GCA" * 30 + "TAA"
        [orf] = [o for o in O.find_orfs(sequence, start_codon_policy="table_starts", min_length=30)
                 if o["start_codon"] == "CTG"]
        expected = "M" + "A" * 30
        self.assertEqual(O.translate_orf_sequence(orf["sequence"]), expected)
        self.assertEqual(P.analyze_translation(orf["sequence"].replace("T", "U"),
                                               initiator_start=True)["protein"], expected)
        self.assertEqual(api.analyze_protein_from_cds(orf["sequence"]).data["protein"], expected)

    def test_initiation_codons_come_from_the_genetic_code(self):
        from analyses import orf as O
        from core.seqtools import get_genetic_code
        self.assertEqual(O.initiation_codons(1), frozenset({"ATG", "CTG", "GTG", "TTG"}))
        for table in self.TABLES:
            with self.subTest(table=table):
                ncbi = get_genetic_code(table).start_codons
                self.assertLessEqual(ncbi, O.initiation_codons(table))
                if table != 1:  # only table 1 adds Synapse's documented GTG policy
                    self.assertEqual(O.initiation_codons(table), ncbi)

    def test_every_initiation_codon_is_decoded_as_methionine(self):
        from analyses import orf as O, protein as P
        for table in self.TABLES:
            for codon in sorted(O.initiation_codons(table)):
                with self.subTest(table=table, codon=codon):
                    cds = codon + "AAATTTTAA"
                    self.assertEqual(O.translate_orf_sequence(cds, genetic_code=table), "MKF")
                    self.assertEqual(P.analyze_translation(cds.replace("T", "U"), initiator_start=True,
                                                           genetic_code=table)["protein"], "MKF")

    def test_internal_codons_keep_their_amino_acid(self):
        from analyses import orf as O, protein as P
        for table, cds, expected in ((1, "ATGCTGGTGTTGTAA", "MLVL"),
                                     (11, "ATGGTGATATTGTAA", "MVIL"),
                                     (2, "ATGATAGTGTAA", "MMV")):
            with self.subTest(table=table):
                self.assertEqual(O.translate_orf_sequence(cds, genetic_code=table), expected)
                self.assertEqual(P.analyze_translation(cds.replace("T", "U"), initiator_start=True,
                                                       genetic_code=table)["protein"], expected)

    def test_a_non_initiation_first_codon_is_not_forced_to_methionine(self):
        from analyses import orf as O, protein as P
        from core import api
        self.assertEqual(O.translate_orf_sequence("GCTAAATTTTAA"), "AKF")
        translation = P.analyze_translation("GCUAAAUUUUAA", initiator_start=True)
        self.assertEqual(translation["protein"], "AKF")
        self.assertFalse(translation["codon_records"][0]["initiator"])
        result = api.analyze_protein_from_cds("GCTAAATTTTAA")
        self.assertEqual(result.data["protein"], "AKF")
        self.assertTrue(any("GCT" in w and "initiation codon" in w for w in result.warnings))

    def test_ambiguous_first_codon_is_an_initiator_only_if_every_reading_is(self):
        from analyses import orf as O, protein as P
        for codon, table, expected in (("WTG", 1, "M"),   # ATG / TTG: both initiation codons
                                       ("YTG", 1, "M"),   # CTG / TTG: both initiation codons
                                       ("NTG", 1, "M"),   # ATG / CTG / GTG / TTG: all initiation codons
                                       ("NNG", 1, "?"),   # includes non-initiation codons, several amino acids
                                       ("RTA", 2, "?")):  # ATA is a table-2 start, GTA (Val) is not
            with self.subTest(codon=codon, table=table):
                self.assertEqual(O.translate_orf_sequence(codon + "AAATAA", genetic_code=table)[0], expected)
                self.assertEqual(P.analyze_translation(codon.replace("T", "U") + "AAAUAA", initiator_start=True,
                                                       genetic_code=table)["protein"][0], expected)

    @unittest.skipUnless(__import__("importlib").util.find_spec("Bio"), "Biopython (reference extra) not installed")
    def test_translation_agrees_with_biopython(self):
        import itertools
        from Bio.Seq import Seq
        from Bio.Data.CodonTable import unambiguous_dna_by_id
        from analyses import orf as O
        codons = ["".join(c) for c in itertools.product("TCAG", repeat=3)]
        for table in self.TABLES:
            ncbi_starts = set(unambiguous_dna_by_id[table].start_codons)
            with self.subTest(table=table, part="internal codons"):
                for codon in codons:
                    reference = str(Seq("ATG" + codon).translate(table=table))[1:]
                    ours = O.translate_orf_sequence("ATG" + codon, genetic_code=table)[1:]
                    self.assertEqual(ours, "" if reference == "*" else reference, codon)
            with self.subTest(table=table, part="first codon"):
                decoded_as_met = {c for c in codons if O.translate_orf_sequence(c + "AAA", genetic_code=table)[:1] == "M"}
                elongation_met = {c for c in codons if str(Seq(c).translate(table=table)) == "M"}
                for codon in ncbi_starts:  # Biopython's cds=True initiator rule
                    self.assertEqual(str(Seq(codon + "AAATAA").translate(table=table, cds=True)), "MK")
                synapse_policy = {"GTG"} if table == 1 else set()
                self.assertEqual(decoded_as_met, ncbi_starts | synapse_policy | elongation_met)


# ----------------------------------------------------------------------
# D2: an ORF's genetic code travels with the ORF into translation
# ----------------------------------------------------------------------
class OrfGeneticCodeHandoffTests(unittest.TestCase):
    MITO = "ATG" + "TGA" * 25 + "TAG"   # TGA = Trp in table 2 (stop in table 1)
    MITO_PROTEIN = "M" + "W" * 25

    def _table2_orf(self):
        from analyses import orf as O
        [orf] = [o for o in O.find_orfs(self.MITO, genetic_code=2, min_length=30, include_partial=False)
                 if o["strand"] == "+" and o["start"] == 0]
        self.assertEqual(orf["genetic_code"], 2)
        return orf

    def test_translation_reads_the_orfs_own_table(self):
        from analyses import orf as O
        orf = self._table2_orf()
        self.assertEqual(O.translate_orf_sequence(orf), self.MITO_PROTEIN)
        self.assertEqual(O.translate_orf_sequence(orf, genetic_code=2), self.MITO_PROTEIN)
        # A bare sequence carries no table: the documented default, table 1, applies.
        self.assertEqual(O.translate_orf_sequence(orf["sequence"]), "M")

    def test_conflicting_table_is_rejected_not_silently_used(self):
        from analyses import orf as O
        from core import api
        orf = self._table2_orf()
        with self.assertRaises(ValueError):
            O.translate_orf_sequence(orf, genetic_code=1)
        with self.assertRaises(ValueError):
            api.analyze_protein_from_cds(orf, genetic_code=1)

    def test_selected_orf_analysis_uses_its_table(self):
        from analyses import orf as O
        self.assertEqual(O.analyze_selected_orf(self._table2_orf())["protein_sequence"], self.MITO_PROTEIN)

    def test_api_protein_from_an_orf_uses_its_table(self):
        from core import api
        result = api.analyze_protein_from_cds(self._table2_orf())
        self.assertEqual(result.data["protein"], self.MITO_PROTEIN)
        self.assertEqual(result.parameters["genetic_code"], 2)
        self.assertEqual(api.analyze_protein_from_cds(self.MITO).data["protein"], "M")  # string: table 1

    def test_api_and_orf_translator_agree_for_every_found_orf(self):
        from analyses import orf as O
        from core import api
        from tests._support import deterministic_dna
        sequence = deterministic_dna(3000, seed=5)
        for table, policy in ((1, "table_starts"), (2, "table_starts"), (4, "atg_and_alternative"), (11, "table_starts")):
            found = api.find_orfs(sequence, genetic_code=table, start_codon_policy=policy, min_length=60,
                                  top_n=200).data["orfs"]
            self.assertTrue(found, table)
            for orf in found:
                with self.subTest(table=table, start=orf["start"], strand=orf["strand"]):
                    protein = O.translate_orf_sequence(orf)
                    self.assertEqual(protein[0], "M")  # every found ORF starts at an initiation codon
                    self.assertEqual(api.analyze_protein_from_cds(orf).data["protein"], protein)


if __name__ == "__main__":
    unittest.main()
