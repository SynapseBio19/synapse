"""Dialog paths behind the first menu: missing and invalid input, cancel, exports.

The rendered-menu tests draw every menu; these drive the variant,
mutation, ORF and literature dialogs to the end of their branches
(Master Prompt 4.2, 11.4, 11.5): empty and invalid answers, operations
before any data exist, files that cannot be read or written, and the
exports into a temporary folder. The literature search is run against a
stub, so no request leaves the machine.
"""

import builtins
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests import _support  # noqa: F401  (blocks network and input())
from tests.test_terminal_ui import _scripted
from ui import style

FIXTURES = Path(__file__).parent / "fixtures"
VCF = FIXTURES / "bcftools_norm" / "norm.split.vcf"


def drive(answers, call):
    """Run a dialog with scripted answers; return (output, result)."""
    output = io.StringIO()
    with mock.patch.object(builtins, "input", _scripted(answers)), contextlib.redirect_stdout(output):
        result = call()
    return output.getvalue(), result


def puc19():
    lines = (FIXTURES / "pUC19_L09137.2.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines if line and not line.startswith(">"))


# ----------------------------------------------------------------------
# Variant analysis
# ----------------------------------------------------------------------
class VariantDialogTests(unittest.TestCase):
    def setUp(self):
        from ui import variant_analysis_menu
        self.menu = variant_analysis_menu

    def test_every_operation_without_variants_says_so(self):
        answers = []
        for choice in ("3", "4", "5", "6", "7", "8", "9", "10", "11"):
            answers += [choice, ""]          # each one pauses after its note
        output, result = drive(answers + ["x", "0"], self.menu.run_variant_analysis_ui)
        self.assertIsNone(result)
        self.assertEqual(output.count("NOTE: No variants loaded."), 9)
        self.assertIn(style.INVALID_SELECTION, output)

    def test_missing_file_empty_paste_and_bad_position(self):
        missing = str(FIXTURES / "does_not_exist.vcf")
        answers = ["1", missing,          # unreadable file
                   "1", "", "",           # paste mode, nothing pasted
                   "2", "chr1", "abc",    # manual variant with a non-integer position
                   "2", "chr1", "",       # ... and with no position
                   "0"]
        output, _ = drive(answers, self.menu.run_variant_analysis_ui)
        self.assertIn("ERROR: Could not read VCF file", output)
        self.assertIn("ERROR: No VCF text supplied.", output)
        # Text is rejected by parse_int, which names the value; an empty answer
        # reaches the dialog's own check.
        self.assertIn("ERROR: Invalid integer value: 'abc'", output)
        self.assertIn("ERROR: Position must be an integer.", output)

    def test_loaded_vcf_through_every_analysis_and_export(self):
        with tempfile.TemporaryDirectory() as folder:
            vcf_out = str(Path(folder) / "out")          # extension added by the dialog
            tsv_out = str(Path(folder) / "variants.tsv")
            answers = [
                "1", str(VCF),
                "5",                                     # statistics
                "8",                                     # cohort / samples
                "4", "500", "", "", "", "",              # filter: QUAL >= 500
                "7", "20\t3\t.\tG\tCT\t999\tPASS\t.", "",  # compare with one pasted record
                "7", "",                                 # compare with nothing pasted
                "3",                                     # normalize without a pipeline sequence
                "9", vcf_out,
                "10", tsv_out,
                "11",
                "12",
                "0",
            ]
            output, result = drive(answers, self.menu.run_variant_analysis_ui)
            written_vcf = Path(vcf_out + ".vcf").read_text()
            written_tsv = Path(tsv_out).read_text()
        self.assertIn("Loaded 29 variant(s).", output)
        self.assertIn("VARIANT STATISTICS", output)
        self.assertIn("COHORT / SAMPLE ANALYSIS", output)
        self.assertIn("XY00001:", output)
        self.assertIn("Passed:", output)
        self.assertIn("VARIANT SET COMPARISON", output)
        self.assertIn("ERROR: No second VCF set supplied.", output)
        self.assertIn("NOTE: No Synapse sequence loaded: variants are trimmed only", output)
        self.assertIn("Normalized 29 variant(s): 0 with a reference sequence", output)
        self.assertIn("OK: VCF exported:", output)
        self.assertIn("OK: TSV exported:", output)
        self.assertIn("VARIANT ANALYSIS SELF-TEST", output)
        self.assertIn("#CHROM\tPOS\tID\tREF\tALT", written_vcf)
        self.assertEqual(len([line for line in written_vcf.splitlines() if not line.startswith("#")]), 29)
        self.assertEqual(len(written_tsv.splitlines()), 30)  # header + 29 rows
        self.assertIsNotNone(result)

    def test_export_into_a_missing_folder_is_an_error(self):
        target = str(FIXTURES / "no_such_folder" / "x.vcf")
        output, _ = drive(["1", str(VCF), "9", target, "10", target[:-3] + "tsv", "0"],
                          self.menu.run_variant_analysis_ui)
        self.assertIn("ERROR: could not write the VCF file", output)
        self.assertIn("ERROR: could not write the TSV file", output)
        self.assertFalse((FIXTURES / "no_such_folder").exists())


# ----------------------------------------------------------------------
# Mutation analysis
# ----------------------------------------------------------------------
REFERENCE = "ATGAAACGTTGGCTGCAAGGTTAA"
VARIANT = "ATGAAATGTTGGCTGCAAGGTTAA"


class MutationDialogTests(unittest.TestCase):
    def setUp(self):
        from ui import mutation_menu
        self.menu = mutation_menu

    def test_statistics_and_report_before_a_comparison(self):
        output, result = drive(["2", "4", "x", "0"], self.menu.run_mutation_analysis_ui)
        self.assertIsNone(result)
        self.assertEqual(output.count("NOTE: Run a sequence comparison first."), 2)
        self.assertIn(style.INVALID_SELECTION, output)

    def test_manual_comparison_statistics_coding_report_and_self_test(self):
        answers = ["1", "Ref", REFERENCE, "Var", VARIANT,   # compare, sequences entered manually
                   "2",
                   "3", "", "",                              # coding: Enter uses the sequences
                   "4",
                   "5",
                   "0"]
        output, result = drive(answers, self.menu.run_mutation_analysis_ui)
        self.assertIn("MUTATION ANALYSIS REPORT", output)
        self.assertIn("c.7C>T", output)
        self.assertIn("MUTATION STATISTICS", output)
        self.assertIn("CODING CONSEQUENCE ANALYSIS", output)
        self.assertIn("Missense:", output)
        self.assertIn("MUTATION SELF-TEST", output)
        self.assertEqual(len(result.mutations), 1)

    def test_cancel_from_the_sequence_choice(self):
        sequences = [("seq1", REFERENCE), ("seq2", VARIANT)]
        output, result = drive(["1", "9", "0", "0"],
                               lambda: self.menu.run_mutation_analysis_ui(sequences=sequences))
        self.assertIn(style.INVALID_SELECTION, output)     # 9 is not listed
        self.assertNotIn("MUTATION ANALYSIS REPORT", output)
        self.assertIsNone(result)

    def test_invalid_sequence_is_reported(self):
        output, _ = drive(["1", "Ref", REFERENCE, "Var", "ATG123", "0"],
                          self.menu.run_mutation_analysis_ui)
        self.assertIn("ERROR:", output)
        self.assertNotIn("MUTATION ANALYSIS REPORT", output)

    def test_frame_warning_uses_the_warning_prefix(self):
        answers = ["3", "Ref", REFERENCE + "A", "Var", VARIANT + "A", "", "", "0"]
        output, _ = drive(answers, self.menu.run_mutation_analysis_ui)
        self.assertIn("WARNING: Coding frame assumption:", output)
        self.assertNotIn("WARNING - coding frame", output)


# ----------------------------------------------------------------------
# ORF analysis
# ----------------------------------------------------------------------
class OrfDialogTests(unittest.TestCase):
    def setUp(self):
        from ui import orf_menu
        self.menu = orf_menu

    def test_every_view_and_the_selection(self):
        answers = ["1", "", "2", "", "3", "", "4", "", "x",
                   "5", "999", "1"]                     # invalid then valid selection
        output, selected = drive(answers, lambda: self.menu.run_orf_analysis_ui(puc19()))
        for heading in ("ORF ANALYSIS", "SELECT ORF", "ORF SELECTED"):
            self.assertIn(heading, output)
        self.assertIn(style.INVALID_SELECTION, output)
        assert selected is not None
        self.assertEqual(selected["rank"], 1)
        self.assertEqual(selected["length"], 861)          # bla, the longest ORF on pUC19
        self.assertEqual(selected["search_parameters"]["min_length_nt"], 75)

    def test_settings_change_every_option_and_reject_bad_values(self):
        answers = ["6",
                   "1", "2",                 # minimum length below 3
                   "1", "abc",               # not a number
                   "1", "78",
                   "2", "3", "4", "5",       # start codons, nested, ranking, partial ORFs
                   "6", "99",                # unknown genetic code
                   "6", "11",
                   "x",
                   "0",
                   "0"]
        output, selected = drive(answers, lambda: self.menu.run_orf_analysis_ui(puc19()))
        self.assertIsNone(selected)
        self.assertEqual(output.count("ERROR: Please enter an integer >= 3."), 2)
        self.assertIn("ERROR: Genetic code not available", output)
        self.assertIn("78 nt here gives ORFfinder's 75 nt", output)
        self.assertIn(">= 78 nt, atg_and_alternative, nested, ranked by heuristic", output)
        self.assertIn("Genetic code:", output)

    def test_restore_defaults_and_legacy(self):
        answers = ["6", "1", "120", "8", "7", "0", "0"]
        output, _ = drive(answers, lambda: self.menu.run_orf_analysis_ui(puc19()))
        self.assertIn(">= 75 nt, atg_only, longest per stop, ranked by length", output)

    def test_selection_on_a_sequence_without_orfs(self):
        output, selected = drive(["5", "0"], lambda: self.menu.run_orf_analysis_ui("ACGT" * 5))
        self.assertIn("NOTE: No ORF candidates available for selection.", output)
        self.assertIsNone(selected)


# ----------------------------------------------------------------------
# Literature search (stubbed provider layer)
# ----------------------------------------------------------------------
class LiteratureDialogTests(unittest.TestCase):
    def setUp(self):
        from analyses import literature as L
        from ui import literature_menu
        self.L, self.menu = L, literature_menu

    def _stub(self, articles, success=True):
        L = self.L

        def search(query):
            run = L.LiteratureSearchRun(query=query, run_id="LIT-TEST")
            run.articles = list(articles)
            run.source_status = {source: "ok" for source in query.sources}
            return L.LiteratureSearchResult(success=success, run=run)
        return search

    def _articles(self):
        L = self.L
        return [
            L.LiteratureArticle(title="Beta-lactamase TEM-1 structure", doi="10.1000/a",
                                authors=["A. Author"], journal="J Test", publication_date="2020-01-01"),
            L.LiteratureArticle(title="Lactase persistence in Europe", pmid="123",
                                authors=["B. Author"], publication_date="2021-06-01"),
        ]

    def test_menu_options_without_a_search(self):
        output, _ = drive(["4", "", "2", "x", "0"], self.menu._interactive_menu)
        self.assertIn("Relevance scores measure query-text matching only", output)
        self.assertIn("SYNAPSE_NCBI_EMAIL", output)
        self.assertIn(style.INVALID_SELECTION, output)

    def test_empty_query_bad_number_and_no_source(self):
        answers = ["1", "",                       # empty query
                   "1", "lactase", "abc",         # maximum results not a number
                   "1", "lactase", "", "9",       # no valid source
                   "0"]
        with mock.patch.object(self.menu, "run_literature_search") as search:
            output, _ = drive(answers, self.menu._interactive_menu)
        search.assert_not_called()
        self.assertIn("ERROR: Query cannot be empty.", output)
        self.assertIn("ERROR: Invalid configuration:", output)
        self.assertIn("ERROR: Select at least one valid source.", output)

    def test_search_detail_and_export(self):
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / "results.json"
            answers = ["1", "lactase", "5", "1,4", "2019", "", "", "date",
                       "x",                         # not a number
                       "7",                         # outside the range
                       "2",                         # detail of result 2
                       "",                          # done with details
                       str(target),
                       "0"]
            with mock.patch.object(self.menu, "run_literature_search", side_effect=self._stub(self._articles())):
                output, _ = drive(answers, self.menu._interactive_menu)
            exported = json.loads(target.read_text(encoding="utf-8"))
        self.assertIn("ERROR: Enter a valid result number.", output)
        self.assertIn("ERROR: Article number is outside the result range.", output)
        self.assertIn("LITERATURE ARTICLE DETAIL", output)
        self.assertIn("OK: Exported:", output)
        titles = json.dumps(exported)
        self.assertIn("Lactase persistence in Europe", titles)
        self.assertIn("Beta-lactamase TEM-1 structure", titles)

    def test_failed_search_and_bad_export_format(self):
        answers = ["1", "lactase", "", "", "", "", "", "",
                   "1", "lactase", "", "", "", "", "", "", "", "out.docx",
                   "0"]
        L = self.L
        failed = L.LiteratureSearchResult(success=False, run=None, message="All sources failed.")
        results = [failed, self._stub(self._articles())(L.LiteratureQuery("lactase"))]
        with mock.patch.object(self.menu, "run_literature_search", side_effect=results):
            output, _ = drive(answers, self.menu._interactive_menu)
        self.assertIn("ERROR: All sources failed.", output)
        self.assertIn("ERROR: Export failed: Supported export formats", output)


if __name__ == "__main__":
    unittest.main()
