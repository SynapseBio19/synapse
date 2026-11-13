"""Execute the features that no test reached before 2026-11-12 (Master Prompt 4.2, T7).

A line-coverage run (coverage 7.10.7, branch mode) showed that the protein
and RNA dialogs ran only as far as their first menu, and that several
export and report functions were never called: the BLAST JSON/CSV/text
exports, the CRISPR CSV/FASTA exports and the base- and prime-editing
reports, and the variant grouping, density and report helpers. Each test
below runs one of them and checks what it produces, not only that it does
not raise.
"""

from __future__ import annotations

import builtins
import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests._support import deterministic_dna
from tests.test_p0_regressions import _blast_xml

PROTEIN = ("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKW"
           "KRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAV")


def _drive(answers, call):
    """Run a dialog with scripted answers; return what it printed."""
    pending = list(answers)
    buffer = io.StringIO()

    def answer(prompt=""):
        buffer.write(prompt)
        if not pending:
            raise AssertionError(f"script exhausted at prompt {prompt!r}")
        return pending.pop(0)

    with mock.patch.object(builtins, "input", answer), contextlib.redirect_stdout(buffer):
        call()
    return buffer.getvalue()


class DialogExecutionTests(unittest.TestCase):
    def test_every_protein_view_is_drawn(self):
        from ui import protein_menu
        output = _drive(["1", "", "2", "", "3", "", "4", "", "5", "", "0"],
                        lambda: protein_menu.run_protein_analysis_from_record(PROTEIN, name="demo"))
        for heading in ("PROTEIN OVERVIEW", "TRANSLATION", "AMINO ACID COMPOSITION",
                        "PHYSICOCHEMICAL", "COMPLETE PROTEIN REPORT"):
            self.assertIn(heading, output.upper())
        self.assertIn("Molecular weight", output)
        self.assertNotIn("Traceback", output)
        self.assertNotIn("ERROR:", output)

    def test_every_rna_view_is_drawn(self):
        from ui import rna_menu
        dna = deterministic_dna(600)
        output = _drive(["", "1", "", "2", "", "3", "", "4", "", "5", "", "0"],
                        lambda: rna_menu.run_rna_analysis_ui(dna, header="demo"))
        upper = output.upper()
        for heading in ("RNA OVERVIEW", "TRANSLATION", "RNA FEATURES", "RNA SEQUENCE"):
            self.assertIn(heading, upper)
        # The sequence view shows the transcript: T is replaced by U.
        self.assertIn(dna[:30].replace("T", "U"), output.replace("\n", "").replace(" ", ""))
        self.assertNotIn("ERROR:", output)


class BlastExportTests(unittest.TestCase):
    def _run(self):
        from analyses import blast as B
        query = B.build_query("A" * 40, name="q")
        return B.parse_blast_xml(_blast_xml(1, 32, 32, 1e-5, "A" * 32, "A" * 32), query=query, rid="T",
                                 parameters=B.BlastParameters(), submitted_at="2026-01-01T00:00:00+00:00")

    def test_json_csv_and_text_exports(self):
        from analyses import blast as B
        run = self._run()
        with tempfile.TemporaryDirectory() as folder:
            data = json.loads(Path(B.export_json(run, Path(folder) / "run.json")).read_text(encoding="utf-8"))
            self.assertEqual(data["metadata"]["rid"], "T")
            self.assertEqual(len(data["hits"]), 1)
            with open(B.export_csv(run, Path(folder) / "run.csv"), encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            raw = Path(B.export_text(run, Path(folder) / "run.txt")).read_bytes()
            # D24: each line ended in CR CR LF before 2026-11-12.
            self.assertNotIn(b"\r\r\n", raw)
            text = raw.decode("utf-8").replace("\r\n", "\n")
        self.assertEqual(text.splitlines()[:2], ["=" * 70, "SYNAPSE BLAST REPORT"])
        self.assertIn("RID:                   T", text)
        # 32 of the 40 query bases are covered by the HSP.
        self.assertIn("Coverage:   80.00%", text)
        self.assertIn("E-value:    1e-05", text)
        self.assertEqual(rows[0]["accession"] if "accession" in rows[0] else rows[0].get("hit_accession"), "X1")


class CrisprExportTests(unittest.TestCase):
    def setUp(self):
        from analyses import crispr as C
        self.C = C
        self.dna = deterministic_dna(500)
        self.result = C.design_guides(self.dna, sequence_name="demo")
        self.assertTrue(self.result.candidates)

    def test_guide_dict_csv_and_fasta(self):
        C = self.C
        first = self.result.candidates[0]
        self.assertEqual(C.guide_to_dict(first)["guide"], first.guide)
        with tempfile.TemporaryDirectory() as folder:
            with open(C.export_csv(self.result, Path(folder) / "guides.csv"), encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), len(self.result.candidates))
            self.assertEqual(rows[0]["guide"], first.guide)
            fasta = Path(C.export_guides_fasta(self.result, Path(folder) / "guides.fa", top_n=3)).read_text(encoding="utf-8")
        records = [line for line in fasta.splitlines() if line.startswith(">")]
        self.assertEqual(len(records), 3)
        self.assertIn(first.guide, fasta)

    def test_mutation_target_and_editing_reports(self):
        C = self.C
        target = C.parse_mutation_target("chr1:120 C>T")
        self.assertEqual((target.chromosome, target.position, target.reference, target.alternative),
                         ("chr1", 120, "C", "T"))
        with self.assertRaises(C.CRISPRAnalysisError):
            C.parse_mutation_target("120 C to T")
        first = self.result.candidates[0]
        self.assertEqual(C.base_editing_window(first, "CBE"), tuple(C.get_base_editor("CBE")["window"]))
        self.assertEqual(C.base_editing_window(first, "CBE", (3, 9)), (3, 9))
        position = self.dna.find("C", 200) + 1
        mutation = C.MutationTarget("demo", position, "C", "T")
        base = C.format_base_editing_report(C.analyze_base_editing(self.result.candidates, mutation))
        self.assertTrue(base.startswith("BASE EDITING ANALYSIS"))
        self.assertIn("No base-editing candidates available.", C.format_base_editing_report([]))
        prime = C.format_prime_editing_report(
            C.design_prime_editing_candidates(self.dna, mutation, self.result.candidates))
        self.assertTrue(prime.startswith("PRIME EDITING CANDIDATES"))
        self.assertIn("No sequence-level prime-editing candidates available.", C.format_prime_editing_report([]))


VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t50\tPASS\t.
chr1\t400\t.\tC\tT\t20\tPASS\t.
chr2\t50\t.\tG\tGA\t.\tPASS\t.
"""


class VariantHelperTests(unittest.TestCase):
    def test_grouping_density_report_and_consequence(self):
        from analyses import variant_analysis as V
        result = V.analyze_vcf_text(VCF, name="demo")
        variants = result.variants
        self.assertEqual(len(variants), 3)
        self.assertEqual(len(V.variant_key_set(variants)), 3)
        self.assertEqual({k: len(v) for k, v in V.variants_by_chromosome(variants).items()}, {"chr1": 2, "chr2": 1})
        self.assertEqual(list(V.variants_by_gene(variants)), ["Unknown"])
        self.assertEqual(list(V.variants_by_consequence(variants)), ["unknown"])
        # Two variants in 1-1000 of chr1: 2 per 1,000 bp.
        self.assertAlmostEqual(V.interval_variant_density(variants, "chr1", 1, 1000), 2.0)
        with self.assertRaises(V.VariantAnalysisError):
            V.interval_variant_density(variants, "chr1", 10, 5)
        self.assertEqual(V.format_variant(variants[0]), "chr1:100 A>G [SNV] QUAL=50.00")
        self.assertTrue(V.format_variant(variants[2]).endswith("QUAL=."))
        report = V.build_variant_report(result)
        self.assertIn("SYNAPSE VARIANT ANALYSIS", report)
        self.assertIn("chr1:400 C>T", report)
        variants[0].consequence = "synonymous_variant&stop_gained"
        self.assertEqual(V.highest_priority_consequence(variants[0]), "stop_gained")
        self.assertIsNone(V.highest_priority_consequence(variants[1]))


if __name__ == "__main__":
    unittest.main()
