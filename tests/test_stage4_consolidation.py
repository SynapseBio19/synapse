"""Stage 4 tests: shared tables, genetic code in protein analysis, core API."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from tests._support import PROJECT_ROOT, deterministic_dna, normalize_for_json


class SharedTableTests(unittest.TestCase):
    def test_complement_tables_are_shared(self):
        from core.seqtools import IUPAC_COMPLEMENT
        from analyses import crispr, mutation, primer
        for module in (crispr, mutation, primer):
            self.assertEqual(module.IUPAC_COMPLEMENT, IUPAC_COMPLEMENT)

    def test_codon_dict_formats(self):
        from core.seqtools import codon_dict
        self.assertEqual(codon_dict(1)["TAA"], "*")
        self.assertEqual(codon_dict(1, rna=True, stop_symbol="STOP")["UGA"], "STOP")
        self.assertEqual(codon_dict(2)["TGA"], "W")
        first = codon_dict(1)
        first["AAA"] = "X"
        self.assertEqual(codon_dict(1)["AAA"], "K")  # new dict per call

    def test_primer_has_single_live_definition_per_name(self):
        import ast
        tree = ast.parse((PROJECT_ROOT / "analyses" / "primer.py").read_text(encoding="utf-8"))
        names = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)]
        duplicates = {n for n in names if names.count(n) > 1}
        # calculate_tm_nearest_neighbor is defined twice on purpose: the v3
        # model is captured under _calculate_tm_nearest_neighbor_v3_raw before
        # the adapter replaces the public name.
        self.assertEqual(duplicates, {"calculate_tm_nearest_neighbor"})


class ProteinGeneticCodeTests(unittest.TestCase):
    def test_translation_with_table_2(self):
        from analyses import protein as P
        rna = "AUG" + "UGA" * 5 + "AGA"
        standard = P.analyze_translation(rna, initiator_start=True)
        mito = P.analyze_translation(rna, initiator_start=True, genetic_code=2)
        self.assertEqual((standard["protein"], standard["genetic_code"]), ("M", 1))
        self.assertEqual((mito["protein"], mito["stop_found"]), ("MWWWWW", True))

    def test_pipeline_uses_orf_genetic_code(self):
        from core.pipeline import SynapsePipeline
        orf_sequence = "ATG" + "TGA" * 30 + "AGA"
        pipeline = SynapsePipeline(orf_sequence, "mito")
        pipeline.set_selected_orf({"sequence": orf_sequence, "rank": 1, "genetic_code": 2})
        result = pipeline.run_protein()
        self.assertEqual(result["protein"], "M" + "W" * 30)


class CoreApiTests(unittest.TestCase):
    SEQUENCE = deterministic_dna(600, seed=11)

    def test_orf_api_matches_module(self):
        from core import api
        from analyses.orf import analyze_orfs_enhanced
        result = api.find_orfs(self.SEQUENCE)
        direct = analyze_orfs_enhanced(self.SEQUENCE)
        self.assertEqual(normalize_for_json(result.data), normalize_for_json(direct))
        self.assertEqual(result.parameters["min_length_nt"], 75)
        self.assertEqual(len(result.input_sequence_sha256), 64)

    def test_restriction_api_uses_record_topology(self):
        from core import api
        from core.models import SequenceRecord
        record = SequenceRecord("p", "ATTC" + "G" * 20 + "GA", topology="circular")
        result = api.restriction_digest(record, ["EcoRI"])
        self.assertEqual(result.parameters["topology"], "circular")

    def test_other_entry_points_run_without_input(self):
        from core import api
        cds = "ATG" + "GCT" * 40 + "TAA"
        self.assertEqual(api.sequence_overview(self.SEQUENCE).module, "sequence")
        self.assertEqual(api.analyze_rna(self.SEQUENCE).module, "rna")
        protein = api.analyze_protein_from_cds(cds)
        self.assertTrue(protein.data["protein"].startswith("MA"))
        comparison = api.compare_sequences(cds, cds[:9] + "C" + cds[10:])  # G10C
        self.assertEqual(comparison.data.statistics.total_mutations, 1)
        guides = api.design_crispr_guides(self.SEQUENCE)
        self.assertEqual(guides.parameters["system"], "SpCas9")
        records = api.load_records(str(PROJECT_ROOT / "examples" / "dna.fasta"))
        self.assertEqual(records[0].identifier, "TEST_8_REFERENCE")

    def test_vcf_api(self):
        from core import api
        vcf = ("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
               "c1\t2\t.\tC\tT\t50\tPASS\t.\nc1\t4\t.\tG\tA\t50\t.\t.\n")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "v.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(vcf)
            plain = api.analyze_vcf(str(path))
            checked = api.analyze_vcf(str(path), contig_references={"c1": "ACGTA"})
        self.assertEqual(plain.data["statistics"].unfiltered, 1)
        self.assertTrue(any("not left-normalized" in w for w in plain.warnings))
        self.assertEqual(checked.data["reference_check"]["mismatched"], 1)  # POS 4 is T, not G
        self.assertEqual(len(checked.parameters["vcf_sha256"]), 64)

    def test_primer_api(self):
        from core import api
        result = api.design_primers(self.SEQUENCE, 100, 400)
        self.assertEqual(result.module, "primer")
        self.assertIn("santalucia_1998", result.limitations[0])


if __name__ == "__main__":
    unittest.main()
