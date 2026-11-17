"""Library functions that FEATURE_STATUS listed as tested but the suite never ran.

Found on 2026-11-17 by checking every FEATURE_STATUS row against line
coverage (Phase 7, package A). Four of these functions have no caller in
the menus (they are part of the module APIs); the three literature
providers are called by the search but were only stubbed out. Each is now
run with checked output. The literature tests replace the HTTP layer
(``_request_json``) with recorded-shape responses; nothing is sent.
"""

import os
import unittest
from types import SimpleNamespace
from unittest import mock
from urllib.parse import parse_qs, urlparse

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import literature as L
from analyses import mutation as M
from analyses import primer as P
from analyses import restriction as R

PUC19_MCS = "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT"  # EcoRI ... HindIII


class DoubleDigestTests(unittest.TestCase):
    def test_two_enzymes_cut_and_fragments_add_up(self):
        sequence = "ATGC" * 10 + PUC19_MCS + "TTAA" * 10
        result = R.run_double_digest(sequence, "EcoRI", "HindIII")
        self.assertEqual((result["status"], result["valid"]), ("OK", True))
        fragments = result["digest"]["fragments"]
        self.assertEqual(len(fragments), 3)
        self.assertEqual(sum(f["length"] for f in fragments), len(sequence))
        self.assertEqual(len(result["gel_bands"]), 3)

    def test_the_same_enzyme_twice_is_refused(self):
        result = R.run_double_digest(PUC19_MCS, "EcoRI", "EcoRI")
        self.assertEqual((result["status"], result["valid"], result["digest"]), ("INVALID", False, None))
        self.assertIn("two different enzymes", result["error"])
        with self.assertRaisesRegex(ValueError, "two different enzymes"):
            R.analyze_double_digest(PUC19_MCS, "BamHI", "BamHI")


def _candidate(tm):
    return SimpleNamespace(thermodynamics=SimpleNamespace(tm_C=tm))


class ThermodynamicFailureDiagnosisTests(unittest.TestCase):
    config = P.ProfessionalPrimerConfig()   # Tm window 57-63 °C

    def test_no_candidates_on_either_side(self):
        reasons = P.diagnose_thermodynamic_failure(forward_candidates=[], reverse_candidates=[],
                                                   pairs=[], config=self.config)
        self.assertEqual(reasons, ["Nearest-neighbor Tm rejected all forward candidates.",
                                   "Nearest-neighbor Tm rejected all reverse candidates."])

    def test_candidates_outside_the_window_and_no_pair(self):
        reasons = P.diagnose_thermodynamic_failure(
            forward_candidates=[_candidate(50.0), _candidate(52.0)],
            reverse_candidates=[_candidate(70.0)], pairs=[], config=self.config)
        self.assertIn("No forward/reverse combination satisfies the thermodynamic pair constraints.", reasons)
        self.assertIn("Forward Tm values are consistently below the configured minimum.", reasons)
        self.assertIn("Reverse Tm values are consistently above the configured maximum.", reasons)

    def test_nothing_to_report_inside_the_window(self):
        reasons = P.diagnose_thermodynamic_failure(
            forward_candidates=[_candidate(60.0)], reverse_candidates=[_candidate(59.0)],
            pairs=[object()], config=self.config)
        self.assertEqual(reasons, [])


class IndelLeftNormalizationTests(unittest.TestCase):
    REFERENCE = "GCTAAAAGC"   # A run at 4-7

    def test_insertion_in_a_run_moves_to_its_left_end(self):
        # An A inserted before position 8 (after the run) is the same event as one before 4.
        self.assertEqual(M.normalize_indel_left(self.REFERENCE, "", 8, "", "A"), (4, "", "A"))

    def test_deletion_in_a_run_moves_to_its_left_end(self):
        self.assertEqual(M.normalize_indel_left(self.REFERENCE, "", 7, "A", ""), (4, "A", ""))

    def test_non_repeated_insertion_stays(self):
        self.assertEqual(M.normalize_indel_left(self.REFERENCE, "", 3, "", "G"), (3, "", "G"))

    def test_substitutions_are_returned_unchanged_and_bad_positions_refused(self):
        self.assertEqual(M.normalize_indel_left(self.REFERENCE, "", 5, "A", "T"), (5, "A", "T"))
        with self.assertRaises(M.MutationAnalysisError):
            M.normalize_indel_left(self.REFERENCE, "", 0, "", "A")


class ProteinComparisonTests(unittest.TestCase):
    def test_mismatches_identity_and_length_change(self):
        result = M.compare_proteins("MKTAYIAK", "MKSAYIA")
        self.assertEqual(result["mismatches"], [{"position": 3, "reference": "T", "alternative": "S"}])
        self.assertEqual((result["reference_length"], result["alternative_length"], result["length_change"]),
                         (8, 7, -1))
        self.assertAlmostEqual(result["identity_percent"], 100.0 * 7 / 8)

    def test_empty_proteins_are_identical(self):
        self.assertEqual(M.compare_proteins("", "")["identity_percent"], 100.0)


class LiteratureProviderTests(unittest.TestCase):
    """The provider functions parse each service's documented response shape."""

    def _run(self, function, response, *args, **kwargs):
        seen = []

        def fake(url, *a, **k):
            seen.append(url)
            return response
        with mock.patch.object(L, "_request_json", side_effect=fake):
            articles = function(*args, **kwargs)
        query = parse_qs(urlparse(seen[0]).query)
        return articles, query

    def test_europe_pmc(self):
        response = {"resultList": {"result": [
            {"title": "Lactase persistence in Europe", "pmid": "123", "doi": "10.1000/x",
             "authorString": "Smith A, Jones B.", "journalTitle": "J Test", "firstPublicationDate": "2021-06-01",
             "abstractText": "An abstract.", "isOpenAccess": "Y", "citedByCount": "7",
             "keywordList": {"keyword": ["lactase", " "]}},
            {"title": "", "pmid": "999"},                        # no title: skipped
        ]}}
        articles, query = self._run(L.search_europe_pmc, response, "lactase", 5,
                                    date_from="2020", date_to="2021-06", open_access_only=True)
        self.assertEqual(len(articles), 1)
        article = articles[0]
        self.assertEqual((article.pmid, article.doi, article.authors), ("123", "10.1000/x", ["Smith A", "Jones B."]))
        self.assertEqual((article.is_open_access, article.citation_count, article.keywords), (True, 7, ["lactase"]))
        self.assertEqual(article.source_urls, ["https://pubmed.ncbi.nlm.nih.gov/123/"])
        search = query["query"][0]
        self.assertIn("(lactase)", search)
        self.assertIn("FIRST_PDATE:[2020-01-01 TO *]", search)
        self.assertIn("FIRST_PDATE:[* TO 2021-06-30]", search)
        self.assertIn("OPEN_ACCESS:Y", search)
        self.assertEqual(query["pageSize"], ["5"])

    def test_openalex(self):
        response = {"results": [
            {"id": "https://openalex.org/W1", "title": "Beta-lactamase TEM-1 structure",
             "doi": "https://doi.org/10.1000/Y", "authorships": [{"author": {"display_name": "A. Author"}}, {}],
             "abstract_inverted_index": {"structure": [1], "The": [0], "resolved": [2]},
             "ids": {"pmid": "https://pubmed.ncbi.nlm.nih.gov/456"},
             "primary_location": {"source": {"display_name": "J Struct"}},
             "publication_date": "2020-01-01", "cited_by_count": 12, "open_access": {"is_oa": False}},
            {"id": "", "title": "no id"},                         # skipped
        ]}
        with mock.patch.dict(os.environ, {"SYNAPSE_OPENALEX_EMAIL": "someone@example.org"}):
            articles, query = self._run(L.search_openalex, response, "tem-1", 3, date_from="2019-05")
        self.assertEqual(len(articles), 1)
        article = articles[0]
        self.assertEqual(article.abstract, "The structure resolved")
        self.assertEqual((article.doi, article.pmid, article.journal), ("10.1000/Y", "456", "J Struct"))
        self.assertEqual((article.authors, article.citation_count, article.is_open_access), (["A. Author"], 12, False))
        self.assertEqual(article.article_id, "openalex:W1")
        self.assertEqual(query["filter"], ["from_publication_date:2019-05-01"])
        self.assertEqual(query["mailto"], ["someone@example.org"])

    def test_crossref(self):
        response = {"message": {"items": [
            {"DOI": "10.1000/z", "title": ["Restriction enzymes revisited"],
             "author": [{"given": "Ada", "family": "Lovelace"}, {"family": "Turing"}, {}],
             "container-title": ["Nucleic Acids Res"], "published-print": {"date-parts": [[2003, 1, 5]]},
             "type": "journal-article", "URL": "https://doi.org/10.1000/z"},
            {"title": []},                                       # no title: skipped
        ]}}
        articles, query = self._run(L.search_crossref, response, "restriction enzyme", 2,
                                    date_from="2000", date_to="2005")
        self.assertEqual(len(articles), 1)
        article = articles[0]
        self.assertEqual((article.doi, article.authors), ("10.1000/z", ["Ada Lovelace", "Turing"]))
        self.assertEqual((article.journal, article.publication_date), ("Nucleic Acids Res", "2003-01-05"))
        self.assertEqual(article.publication_types, ["journal-article"])
        self.assertEqual(query["filter"], ["from-pub-date:2000,until-pub-date:2005"])
        self.assertEqual(query["rows"], ["2"])


if __name__ == "__main__":
    unittest.main()
