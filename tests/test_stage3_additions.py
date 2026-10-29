"""Stage 3 tests: P1 scientific improvements and reproducibility features."""

from __future__ import annotations

import gzip
import io
import json
import math
import tempfile
import unittest
import urllib.error
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest import mock

from tests import _support  # noqa: F401
from tests._support import reverse_complement


# ----------------------------------------------------------------------
# ORF: genetic code option (P1-9)
# ----------------------------------------------------------------------
class OrfGeneticCodeTests(unittest.TestCase):
    MITO = "ATG" + "TGA" * 30 + "AGA"  # TGA = Trp and AGA = stop in table 2

    def test_table1_unchanged_and_table2_reads_through_tga(self):
        from analyses import orf as O
        self.assertEqual(O.analyze_orfs(self.MITO, include_partial=False)["orfs"], [])
        result = O.analyze_orfs(self.MITO, include_partial=False, genetic_code=2)
        self.assertEqual([(o["start"], o["stop"], o["genetic_code"]) for o in result["orfs"]], [(0, 95, 2)])
        self.assertEqual(result["parameters"]["genetic_code_table"], 2)
        self.assertEqual(result["parameters"]["stop_codons"], ["AGA", "AGG", "TAA", "TAG"])

    def test_table1_parameters_keep_previous_description(self):
        from analyses import orf as O
        params = O.analyze_orfs("ATG" + "AAA" * 30 + "TAA")["parameters"]
        self.assertEqual(params["genetic_code"], "NCBI translation table 1 (standard)")
        self.assertEqual(params["stop_codons"], ["TAA", "TAG", "TGA"])

    def test_translation_uses_selected_table(self):
        from analyses import orf as O
        self.assertEqual(O.translate_orf_sequence("ATGTGATAA"), "M")
        self.assertEqual(O.translate_orf_sequence("ATGTGATAA", genetic_code=2), "MW")
        # ATA is a start codon in table 2 and is decoded as initiator Met.
        self.assertEqual(O.translate_orf_sequence("ATATGATAA", genetic_code=2), "MW")
        selected = O.analyze_selected_orf(O.analyze_orfs_enhanced(self.MITO, genetic_code=2)["orfs"][0])
        self.assertEqual(selected["protein_sequence"], "M" + "W" * 30)

    def test_start_policies_for_tables(self):
        from analyses import orf as O
        self.assertEqual(O.start_codons_for(2, "atg_only"), {"ATG"})
        self.assertEqual(O.start_codons_for(1, "atg_and_alternative"), {"ATG", "GTG", "TTG"})
        self.assertEqual(O.start_codons_for(1, "table_starts"), {"ATG", "CTG", "TTG"})
        self.assertIn("ATA", O.start_codons_for(2, "atg_and_alternative"))
        with self.assertRaises(ValueError):
            O.start_codons_for(1, "nonsense")
        with self.assertRaises(ValueError):
            O.analyze_orfs("ATGTAA", genetic_code=99)


# ----------------------------------------------------------------------
# Mutation: coding-frame assumptions (P1-10)
# ----------------------------------------------------------------------
class MutationFrameTests(unittest.TestCase):
    def test_plausible_cds_has_no_warnings(self):
        from analyses import mutation as M
        self.assertEqual(M.assess_coding_frame("ATGAAACCCTAA"), [])

    def test_non_coding_frame_is_flagged(self):
        from analyses import mutation as M
        warnings = M.assess_coding_frame("CCCTAAGGGA")
        self.assertEqual(len(warnings), 4)
        self.assertTrue(any("in-frame stop" in w for w in warnings))

    def test_orf_reference_cds_minus_strand(self):
        from analyses import mutation as M
        cds = "ATGAAACCCTAA"
        reference = "GG" + reverse_complement(cds) + "TT"
        variant = reference[:5] + ("A" if reference[5] != "A" else "C") + reference[6:]
        orf = {"sequence": cds, "rank": 1, "strand": "-", "frame": -1,
               "start": 2, "stop": 13, "length": 12}
        ref_cds, var_cds, text = M.orf_reference_cds(orf, reference, variant)
        self.assertEqual(ref_cds, cds)
        self.assertEqual(var_cds, reverse_complement(variant[2:14]))
        self.assertIn("- strand", text)
        _, var_none, _ = M.orf_reference_cds(orf, reference, reference + "A")
        self.assertIsNone(var_none)

    def test_long_alignment_note(self):
        from analyses import mutation as M
        short = M.compare_sequences("ACGT" * 10, "ACGA" + "ACGT" * 9)
        self.assertFalse(any("aligned globally" in n for n in short.notes))
        reference = "ACGT" * 520
        long_result = M.compare_sequences(reference, reference[:100] + "T" + reference[101:])
        self.assertTrue(any("aligned globally" in n for n in long_result.notes))


# ----------------------------------------------------------------------
# Reproducible analysis report (P1-12)
# ----------------------------------------------------------------------
class AnalysisReportTests(unittest.TestCase):
    def _pipeline(self):
        from core.pipeline import SynapsePipeline
        pipeline = SynapsePipeline("ATGAAATTTTAA", "demo")
        pipeline.store("crispr", {"guides": {"b", "a"}, "score": math.nan},
                       parameters={"system": "SpCas9"}, warnings=["w1"], limitations=["l1"])
        return pipeline

    def test_report_contents_are_json_safe_and_hashed(self):
        from core.report import build_analysis_report, canonical_json
        report = build_analysis_report(self._pipeline(), generated_at="2026-01-01T00:00:00+00:00")
        json.dumps(report, allow_nan=False)
        entry = report["analyses"][0]
        self.assertEqual(entry["data"], {"guides": ["a", "b"], "score": "NaN"})
        self.assertEqual(entry["parameters"], {"system": "SpCas9"})
        self.assertEqual(entry["input_sequence_sha256"], report["input"]["sequence_sha256"])
        import hashlib
        self.assertEqual(entry["data_sha256"],
                         hashlib.sha256(canonical_json(entry["data"]).encode()).hexdigest())

    def test_export_json_and_markdown(self):
        from core.report import export_analysis_report
        pipeline = self._pipeline()
        with tempfile.TemporaryDirectory() as tmp:
            json_path = export_analysis_report(pipeline, str(Path(tmp) / "r.json"))
            md_path = export_analysis_report(pipeline, str(Path(tmp) / "sub" / "r.md"))
            data = json.loads(json_path.read_text(encoding="utf-8"))
            markdown = md_path.read_text(encoding="utf-8")
            with self.assertRaises(ValueError):
                export_analysis_report(pipeline, str(Path(tmp) / "r.txt"))
        self.assertEqual(data["analysis_count"], 1)
        self.assertIn("### crispr #1", markdown)
        self.assertIn("not been experimentally validated", markdown)


# ----------------------------------------------------------------------
# Restriction: circular topology
# ----------------------------------------------------------------------
class CircularRestrictionTests(unittest.TestCase):
    def test_site_across_origin(self):
        from analyses import restriction as R
        sequence = "ATTC" + "G" * 20 + "GA"   # GAATTC spans the end->start junction
        self.assertEqual(R.find_restriction_site_matches(sequence, "EcoRI"), [])
        sites = R.find_restriction_site_matches(sequence, "EcoRI", topology="circular")
        self.assertEqual(len(sites), 1)
        self.assertTrue(sites[0]["wraps_origin"])
        self.assertEqual(sites[0]["position"], len(sequence) - 2)
        cuts = R.get_cut_positions(sequence, ["EcoRI"], topology="circular")
        self.assertEqual(len(cuts), 1)


# ----------------------------------------------------------------------
# VCF: header cardinality, PL remapping, gzip, REF check
# ----------------------------------------------------------------------
VCF = """##fileformat=VCFv4.3
##INFO=<ID=XA,Number=A,Type=Integer,Description="per ALT">
##INFO=<ID=XR,Number=R,Type=Integer,Description="per allele, incl. REF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t3\t.\tG\tA,T\t50\t.\tXA=7,8;XR=1,2,3\tGT:AD:PL\t1/2:10,5,6:90,80,70,60,50,40
"""


class VcfTests(unittest.TestCase):
    def test_header_numbers(self):
        from analyses import variant_analysis as V
        numbers = V.parse_vcf_header_numbers(VCF.splitlines())
        self.assertEqual(numbers["INFO"], {"XA": "A", "XR": "R"})
        self.assertEqual(numbers["FORMAT"]["PL"], "G")

    def test_split_uses_header_cardinality(self):
        from analyses import variant_analysis as V
        variants, _, samples = V.parse_vcf(VCF)
        self.assertEqual(samples, ["S1"])
        self.assertEqual([v.alternative for v in variants], ["A", "T"])
        first, second = variants
        self.assertEqual(str(first.info["XA"]), "7")
        self.assertEqual([str(x) for x in second.info["XR"]], ["1", "3"])
        # PL order for alleles 0,1,2: 0/0,0/1,1/1,0/2,1/2,2/2 -> F(j,k)=k(k+1)/2+j
        self.assertEqual(first.sample_data["S1"]["PL"], "90,80,70")
        self.assertEqual(second.sample_data["S1"]["PL"], "90,60,40")
        self.assertEqual(second.sample_data["S1"]["AD"], "10,6")
        self.assertEqual(first.sample_data["S1"]["GT"], "1/0")

    def test_gzip_file_and_unfiltered_statistics(self):
        from analyses import variant_analysis as V
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(VCF)
            variants, _, _ = V.parse_vcf_file(str(path))
        self.assertEqual(len(variants), 2)
        stats = V.calculate_variant_statistics(variants)
        self.assertEqual(stats.unfiltered, 2)
        self.assertEqual(stats.failed_filters, 0)

    def test_reference_allele_check(self):
        from analyses import variant_analysis as V
        variants, _, _ = V.parse_vcf(VCF)
        ok = V.check_reference_alleles(variants, {"chr1": "AAGAA"})
        self.assertEqual((ok["checked"], ok["mismatched"]), (2, 0))
        bad = V.check_reference_alleles(variants, {"chr1": "AACAA"})
        self.assertEqual(bad["mismatched"], 2)
        self.assertIn("REF_MISMATCH", variants[0].annotation["validation_flags"])
        unmapped = V.check_reference_alleles(variants, {"chr2": "A"})
        self.assertEqual(unmapped["not_checked"], 2)


# ----------------------------------------------------------------------
# Literature: date bounds, retry, sequence context
# ----------------------------------------------------------------------
class LiteratureTests(unittest.TestCase):
    def test_date_bounds(self):
        from analyses import literature as L
        self.assertEqual(L._date_bound("2020", False), "2020-01-01")
        self.assertEqual(L._date_bound("2020", True), "2020-12-31")
        self.assertEqual(L._date_bound("2020-02", True), "2020-02-29")
        self.assertEqual(L._date_bound("2021-02", True), "2021-02-28")
        self.assertEqual(L._date_bound("2021-03-04", True), "2021-03-04")
        self.assertIsNone(L._date_bound(None, True))

    def test_retry_on_503_then_success(self):
        from analyses import literature as L

        class FakeResponse(io.BytesIO):
            headers = mock.Mock(get_content_charset=lambda: "utf-8")

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

        error = urllib.error.HTTPError("https://example.org/x", 503, "busy", {"Retry-After": "1"}, None)
        calls = [error, FakeResponse(b"ok")]

        def fake_urlopen(request, timeout=None):
            item = calls.pop(0)
            if isinstance(item, Exception):
                raise item
            return item

        sleeps = []
        with mock.patch.object(L.urllib.request, "urlopen", fake_urlopen), \
                mock.patch.object(L.time, "sleep", sleeps.append):
            self.assertEqual(L._request_text("https://example.org/x"), "ok")
        self.assertIn(1.0, sleeps)

    def test_non_retryable_error_is_reported(self):
        from analyses import literature as L
        error = urllib.error.HTTPError("https://example.org/y", 404, "nf", {}, None)
        with mock.patch.object(L.urllib.request, "urlopen", side_effect=error), \
                mock.patch.object(L.time, "sleep", lambda s: None):
            with self.assertRaises(L.LiteratureSourceError):
                L._request_text("https://example.org/y")

    def test_sequence_context_suggestion(self):
        from analyses import literature as L
        context = L.build_sequence_literature_context("BRCA1|Homo sapiens topology=circular 1234")
        self.assertEqual(context["suggested_query"], "BRCA1 Homo sapiens")
        self.assertEqual(L.build_sequence_literature_context(""), {})


# ----------------------------------------------------------------------
# BLAST: reported search parameters and database date
# ----------------------------------------------------------------------
class BlastMetadataTests(unittest.TestCase):
    def test_reported_parameters(self):
        from analyses import blast as B
        root = ET.fromstring(
            "<BlastOutput><BlastOutput_param><Parameters>"
            "<Parameters_matrix>BLOSUM62</Parameters_matrix>"
            "<Parameters_expect>10</Parameters_expect>"
            "</Parameters></BlastOutput_param></BlastOutput>"
        )
        self.assertEqual(B._parse_parameters(root), {"matrix": "BLOSUM62", "expect": "10"})

    def test_blastdbcmd_info_date(self):
        from analyses import blast as B
        info = B._parse_blastdbcmd_info(
            "Database: synthetic test\n"
            "\t12 sequences; 3,456 total bases\n\n"
            "Date: Oct 1, 2026  10:00 AM\tLongest sequence: 900 bases\n\n"
            "BLASTDB Version: 5\n"
        )
        self.assertEqual(info["database_date"], "Oct 1, 2026  10:00 AM")
        self.assertEqual((info["sequence_count"], info["total_length"]), (12, 3456))
        self.assertEqual(info["blastdb_version"], "5")  # existing contract: string


if __name__ == "__main__":
    unittest.main()
