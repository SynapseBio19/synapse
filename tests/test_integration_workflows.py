"""Cross-module workflows (Master Prompt 10.1-10.4, 10.7, 10.8).

Each test runs a complete workflow through the non-interactive layers
(``core.pipeline`` and ``core.api``) and checks what 10.8 asks for:
selections are transferred, coordinates stay consistent, parameters are
preserved, results are stored, invalidated results disappear, and the
report reflects the analysis that was actually run.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import logging
import os
import tempfile
import unittest

from tests import _support  # noqa: F401
from tests._support import deterministic_dna


def _revcomp(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.CRITICAL)
        self._stdout = contextlib.redirect_stdout(io.StringIO())
        self._stdout.__enter__()

    def tearDown(self):
        self._stdout.__exit__(None, None, None)
        logging.disable(logging.NOTSET)


# ----------------------------------------------------------------------
# 10.1 DNA -> ORF -> selection -> RNA -> protein -> storage -> report
# ----------------------------------------------------------------------
class DnaToProteinWorkflowTests(_Quiet):
    CDS = "ATG" + "GCT" * 10 + "TGG" + "AAA" * 10 + "TAA"

    def _pipeline(self, cds, genetic_code=1):
        from core.pipeline import SynapsePipeline
        dna = "GC" * 15 + _revcomp(cds) + "CG" * 15
        pipeline = SynapsePipeline(dna, "workflow")
        analysis = pipeline.run_orf(min_length=60, genetic_code=genetic_code)
        minus = [orf for orf in analysis["orfs"] if orf["strand"] == "-"]
        self.assertEqual(len(minus), 1)
        return pipeline, dna, minus[0]

    def test_minus_strand_orf_reaches_rna_and_protein_with_its_coordinates(self):
        pipeline, dna, orf = self._pipeline(self.CDS)
        # Forward-strand, 0-based inclusive coordinates point at the CDS.
        self.assertEqual(_revcomp(dna[orf["start"]:orf["stop"] + 1]), self.CDS)
        self.assertEqual(orf["sequence"], self.CDS)
        pipeline.set_selected_orf(orf)
        pipeline.run_rna()
        pipeline.run_protein()
        context = pipeline.context
        self.assertEqual(context.rna_sequence, self.CDS.replace("T", "U"))
        self.assertEqual(context.protein_sequence, "M" + "A" * 10 + "W" + "K" * 10)
        self.assertEqual((context.rna_source, context.protein_source), ("Selected ORF", "Selected ORF"))
        self.assertEqual({m: len(r) for m, r in context.results.items()},
                         {"orf": 1, "rna": 1, "protein": 1})

    def test_the_orf_genetic_code_reaches_the_protein(self):
        # TGA is a stop in table 1 and tryptophan in table 2.
        cds = "ATG" + "GCT" * 10 + "TGA" + "AAA" * 10 + "TAA"
        pipeline, _, orf = self._pipeline(cds, genetic_code=2)
        self.assertEqual(orf["genetic_code"], 2)
        pipeline.set_selected_orf(orf)
        pipeline.run_protein()
        self.assertEqual(pipeline.context.protein_sequence, "M" + "A" * 10 + "W" + "K" * 10)

    def test_selecting_another_orf_removes_stale_results_everywhere(self):
        from core.report import build_analysis_report
        pipeline, _, orf = self._pipeline(self.CDS)
        pipeline.set_selected_orf(orf)
        pipeline.run_rna()
        pipeline.run_protein()
        other = dict(orf, rank=99)
        pipeline.set_selected_orf(other)
        context = pipeline.context
        self.assertIsNone(pipeline.context.latest("rna"))
        self.assertIsNone(pipeline.context.latest("protein"))
        self.assertIsNone(context.rna_sequence)
        self.assertIsNone(context.protein_sequence)
        modules = {entry["module"] for entry in build_analysis_report(pipeline)["analyses"]}
        self.assertEqual(modules, {"orf"})


# ----------------------------------------------------------------------
# 10.2 DNA -> restriction -> primer design -> storage -> report
# ----------------------------------------------------------------------
class RestrictionPrimerWorkflowTests(_Quiet):
    def test_restriction_and_primer_results_share_one_coordinate_frame(self):
        from analyses.restriction import analyze_restriction
        from core import api
        from core.pipeline import SynapsePipeline
        from core.report import build_analysis_report
        base = deterministic_dna(500)
        dna = base[:240] + "GAATTC" + base[246:]
        pipeline = SynapsePipeline(dna, "workflow")
        restriction = analyze_restriction(dna, ["EcoRI"], "linear")
        pipeline.context.store_restriction_analysis(restriction)
        cut = restriction["digest"]["cut_positions"]
        self.assertEqual(cut, [241])  # G^AATTC: one base before the cut
        design = api.design_primers(dna, 150, 350)
        pairs = [getattr(item, "pair", item) for item in design.data.ranked_pairs]
        self.assertTrue(pairs)
        for pair in pairs:
            self.assertEqual(dna[pair.forward.start:pair.forward.end], pair.forward.sequence)
            self.assertEqual(_revcomp(dna[pair.reverse.start:pair.reverse.end]), pair.reverse.sequence)
            self.assertEqual(dna[pair.amplicon_start:pair.amplicon_end], pair.amplicon_sequence)
            self.assertTrue(150 <= pair.amplicon_start and pair.amplicon_end <= 350)
        pipeline.context.store_primer_analysis(design.data)
        report = build_analysis_report(pipeline)
        entries = {entry["module"]: entry for entry in report["analyses"]}
        self.assertEqual(set(entries), {"restriction", "primer"})
        self.assertEqual(entries["restriction"]["parameters"]["enzymes"], ["EcoRI"])
        self.assertIn("enzyme_database_version", entries["restriction"]["parameters"])
        self.assertEqual((entries["primer"]["parameters"]["target_start"],
                          entries["primer"]["parameters"]["target_end"]), (150, 350))
        self.assertEqual(entries["primer"]["warnings"], list(design.data.warnings))


# ----------------------------------------------------------------------
# 10.3 DNA -> ORF selection -> CRISPR design restricted to the ORF
# ----------------------------------------------------------------------
class CrisprWorkflowTests(_Quiet):
    def test_guides_for_a_selected_orf_lie_in_it_and_match_the_dna(self):
        from core import api
        dna = deterministic_dna(900)
        orf = api.find_orfs(dna).data["orfs"][0]
        result = api.design_crispr_guides(dna, selected_orf=orf)
        design = result.data
        self.assertEqual((design.region.start, design.region.end), (orf["start"] + 1, orf["stop"] + 1))
        self.assertTrue(design.candidates)
        for guide in design.candidates:
            target = dna[guide.guide_start - 1:guide.guide_end]  # 1-based inclusive
            self.assertEqual(target if guide.strand == "+" else _revcomp(target), guide.guide)
            self.assertTrue(design.region.start <= guide.guide_start <= guide.guide_end <= design.region.end)
            self.assertTrue(guide.guide_start <= guide.cut_site <= guide.guide_end)
        self.assertTrue(result.parameters["orf_restricted"])
        self.assertIn("not a genome", " ".join(result.limitations))


# ----------------------------------------------------------------------
# 10.4 VCF -> parsing -> reference check -> normalization -> statistics
# ----------------------------------------------------------------------
class VariantWorkflowTests(_Quiet):
    def test_reference_mismatch_is_reported_and_provenance_kept(self):
        from core import api
        reference = deterministic_dna(200)
        good_alt = "A" if reference[49] != "A" else "C"
        wrong_ref = "T" if reference[99] != "T" else "G"
        alt = "C" if wrong_ref != "C" else "A"
        text = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                f"chr1\t50\tmatch\t{reference[49]}\t{good_alt}\t50\tPASS\t.\n"
                f"chr1\t100\tmismatch\t{wrong_ref}\t{alt}\t50\tPASS\t.\n")
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "workflow.vcf")
            with open(path, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
            result = api.analyze_vcf(path, contig_references={"chr1": reference})
            digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        self.assertEqual(result.parameters["vcf_sha256"], digest)
        self.assertIn("1 REF allele(s) do not match the reference.", result.warnings)
        self.assertEqual(len(result.data["variants"]), 2)


# ----------------------------------------------------------------------
# 10.7 analysis -> storage -> parameters -> hashes -> JSON/Markdown export
# ----------------------------------------------------------------------
class ReproducibleReportWorkflowTests(_Quiet):
    def test_exported_report_is_traceable_to_input_and_data(self):
        from core.pipeline import SynapsePipeline
        from core.report import canonical_json, export_analysis_report
        dna = deterministic_dna(600)
        pipeline = SynapsePipeline(dna, "report-workflow")
        pipeline.run_orf()
        pipeline.run_restriction(["EcoRI", "BamHI", "AluI"])
        with tempfile.TemporaryDirectory() as directory:
            json_path = export_analysis_report(pipeline, os.path.join(directory, "r.json"))
            md_path = export_analysis_report(pipeline, os.path.join(directory, "r.md"))
            report = json.loads(json_path.read_text(encoding="utf-8"))
            markdown = md_path.read_text(encoding="utf-8")
        self.assertEqual(report["input"]["sequence_sha256"], hashlib.sha256(dna.encode()).hexdigest())
        self.assertEqual(report["analysis_count"], 2)
        for entry in report["analyses"]:
            recomputed = hashlib.sha256(canonical_json(entry["data"]).encode("utf-8")).hexdigest()
            self.assertEqual(entry["data_sha256"], recomputed)
            self.assertEqual(entry["input_sequence_sha256"], report["input"]["sequence_sha256"])
            self.assertTrue(entry["parameters"], entry["module"])
            self.assertIn(f"### {entry['module']} #1", markdown)
        self.assertIn(report["software"]["synapse_version"], markdown)


if __name__ == "__main__":
    unittest.main()
