"""Tests for the Stage 2 foundation: genetic codes, models, file loading, pipeline store."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tests._support import FIXTURES_DIR, PROJECT_ROOT
from ui import mutation_menu, primer_menu, variant_analysis_menu


class GeneticCodeTests(unittest.TestCase):
    def test_table1_matches_independent_reference_and_module_tables(self):
        # Until Stage 3 four independently written tables existed; table 1 was
        # verified against all of them. Since Stage 4 the modules derive their
        # tables from core.seqtools, so the cross-check uses a preserved copy
        # of one of the hand-written tables (tests/fixtures).
        import sys
        from core.seqtools import get_genetic_code
        from analyses import orf, mutation, rna
        from data.genetic_code import CODON_TABLE
        sys.path.insert(0, str(FIXTURES_DIR))
        from standard_code_reference import STANDARD_CODE

        reference = get_genetic_code(1).codons
        self.assertEqual(reference, {k: ("*" if v == "STOP" else v) for k, v in STANDARD_CODE.items()})

        def normalized(table):
            return {k.replace("U", "T"): ("*" if v in ("STOP", "*") else v) for k, v in table.items()}

        for table in (orf.GENETIC_CODE, mutation.CODON_TABLE, rna.RNA_CODON_TABLE, CODON_TABLE):
            self.assertEqual(reference, normalized(table))

    def test_known_differences_of_other_tables(self):
        from core.seqtools import get_genetic_code
        vertebrate_mito = get_genetic_code(2)
        self.assertEqual(vertebrate_mito.codons["AGA"], "*")
        self.assertEqual(vertebrate_mito.codons["TGA"], "W")
        self.assertEqual(vertebrate_mito.codons["ATA"], "M")
        self.assertEqual(get_genetic_code(4).codons["TGA"], "W")
        self.assertIn("GTG", get_genetic_code(11).start_codons)
        self.assertNotIn("GTG", get_genetic_code(1).start_codons)

    def test_translation_rules(self):
        from core.seqtools import get_genetic_code
        code = get_genetic_code(11)
        self.assertEqual(code.translate("GTGAAATTTTAA"), "VKF")
        self.assertEqual(code.translate("GTGAAATTTTAA", initiator_start=True), "MKF")
        self.assertEqual(code.translate_codon("GCN"), "A")
        self.assertEqual(code.translate_codon("NNN"), "X")

    def test_unknown_table_rejected(self):
        from core.seqtools import get_genetic_code
        with self.assertRaises(ValueError):
            get_genetic_code(99)


class IntervalTests(unittest.TestCase):
    def test_conversions_and_extract(self):
        from core.models import Interval
        parent = "AAACCCGGGTTT"
        interval = Interval.from_1based(4, 6, "-")
        self.assertEqual((interval.start, interval.end, interval.length), (3, 6, 3))
        self.assertEqual(interval.to_1based(), (4, 6))
        self.assertEqual(interval.extract(parent), "GGG")
        self.assertEqual(Interval.from_0based_inclusive(3, 5, "+").extract(parent), "CCC")

    def test_invalid_interval(self):
        from core.models import Interval
        with self.assertRaises(ValueError):
            Interval(5, 2)


class LoaderTests(unittest.TestCase):
    def test_fasta_provenance_and_topology_token(self):
        from fileio.loader import load_sequences
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "p.fasta"
            path.write_text(">pX topology=circular\nACGTACGT\n>lin\nGGGG\n", encoding="utf-8")
            records = load_sequences(str(path))
        self.assertEqual([r["topology"] for r in records], ["circular", "linear"])
        self.assertEqual(records[0]["source_format"], "fasta")
        self.assertEqual(len(records[0]["source_sha256"]), 64)

    def test_genbank_fixture(self):
        from fileio.loader import load_sequences
        records = load_sequences(str(FIXTURES_DIR / "synthetic_circular.gb"))
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record["topology"], "circular")
        self.assertEqual(len(record["sequence"]), 120)
        self.assertEqual(record["notes"], [])
        cds = [f for f in record["features"] if f["key"] == "CDS"][0]
        location = cds["location_parsed"]
        self.assertEqual(location["strand"], "-")
        self.assertEqual(location["spans_1based"], [(10, 30), (40, 60)])
        self.assertEqual(cds["qualifiers"]["product"], "synthetic test protein")
        gene = [f for f in record["features"] if f["key"] == "gene"][0]
        self.assertTrue(gene["location_parsed"]["partial"])

    def test_example_file_still_loads(self):
        from fileio.loader import load_sequences
        records = load_sequences(str(PROJECT_ROOT / "examples" / "dna.fasta"))
        self.assertEqual(records[0]["header"], "TEST_8_REFERENCE")

    def test_unknown_format_rejected(self):
        from fileio.loader import load_sequences
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.txt"
            path.write_text("hello\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                load_sequences(str(path))


class PipelineStoreTests(unittest.TestCase):
    def test_store_and_invalidation(self):
        from core.pipeline import SynapsePipeline
        pipeline = SynapsePipeline("ATGAAATTTTAA", "test")
        pipeline.store("crispr", {"x": 1}, parameters={"system": "SpCas9"})
        pipeline.store("mutation", {"y": 2})
        result = pipeline.latest("crispr")
        self.assertEqual(result.parameters["system"], "SpCas9")
        self.assertEqual(len(result.input_sequence_sha256), 64)
        pipeline.set_selected_orf({"sequence": "ATGAAATTTTAA", "rank": 1})
        self.assertIsNone(pipeline.latest("crispr"))      # ORF-dependent: invalidated
        self.assertIsNotNone(pipeline.latest("mutation"))  # independent: kept

    def test_primer_store_path(self):
        from core.pipeline import SynapsePipeline
        pipeline = SynapsePipeline("ACGT" * 50, "test")
        primer_menu._primer_ui_store_result(pipeline, object())
        self.assertIsNotNone(pipeline.latest("primer"))

    def test_mutation_and_variant_store_paths(self):
        from core.pipeline import SynapsePipeline
        pipeline = SynapsePipeline("ACGT" * 50, "test")
        mutation_menu._store_pipeline_result(pipeline, "mutation", {"a": 1})
        variant_analysis_menu._store_variant_pipeline_result(pipeline, "variant", {"b": 2})
        self.assertEqual(pipeline.latest("mutation").data, {"a": 1})
        self.assertEqual(pipeline.latest("variant").data, {"b": 2})

    def test_attach_record_requires_matching_sequence(self):
        from core.pipeline import SynapsePipeline
        from core.models import SequenceRecord
        pipeline = SynapsePipeline("ACGT", "t")
        with self.assertRaises(ValueError):
            pipeline.attach_record(SequenceRecord("x", "GGGG"))
        pipeline.attach_record(SequenceRecord("x", "ACGT", topology="circular"))
        self.assertEqual(pipeline.context.topology, "circular")


if __name__ == "__main__":
    unittest.main()
