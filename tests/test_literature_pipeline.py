"""PubMed search, the combined search, metadata enrichment and the exports.

Phase 7, package D. These paths make network requests in normal use; here
the HTTP layer (`_request_json`, `_request_text`) or the provider functions
are replaced by stand-ins with recorded-shape responses, so the parsing,
merging, partial-failure handling, enrichment and export logic run without
the network. The records are synthetic.
"""

import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock
from urllib.parse import parse_qs, urlparse

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import literature as L
# Since Phase 10 the implementation lives in analyses._literature; a stand-in
# must replace a helper in the module that calls it.
from analyses._literature import enrichment, providers, search

PUBMED_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
 <PubmedArticle>
  <MedlineCitation><PMID>111</PMID>
   <Article>
    <Journal><Title>Journal of Tests</Title><JournalIssue><PubDate><Year>2020</Year><Month>Mar</Month><Day>04</Day></PubDate></JournalIssue></Journal>
    <ArticleTitle>Lactase <i>persistence</i> in Europe</ArticleTitle>
    <Abstract><AbstractText Label="BACKGROUND">Why.</AbstractText><AbstractText Label="UNASSIGNED">How.</AbstractText></Abstract>
    <AuthorList><Author><LastName>Smith</LastName><ForeName>Ann</ForeName></Author><Author><CollectiveName>Test Consortium</CollectiveName></Author></AuthorList>
    <PublicationTypeList><PublicationType>Journal Article</PublicationType></PublicationTypeList>
   </Article>
  </MedlineCitation>
  <PubmedData><ArticleIdList><ArticleId IdType="doi">10.1000/abc</ArticleId><ArticleId IdType="pmc">PMC9</ArticleId></ArticleIdList></PubmedData>
 </PubmedArticle>
 <PubmedArticle><MedlineCitation><PMID></PMID></MedlineCitation></PubmedArticle>
</PubmedArticleSet>"""


class PubMedSearchTests(unittest.TestCase):
    def test_search_parses_titles_authors_abstract_ids_and_dates(self):
        urls = []

        def fake_json(url, *a, **k):
            urls.append(url)
            return {"esearchresult": {"idlist": ["111", "222"]}}

        with mock.patch.object(providers, "_request_json", side_effect=fake_json), \
                mock.patch.object(providers, "_request_text", return_value=PUBMED_XML):
            articles = L.search_pubmed("lactase", 5, date_from="2019-01", date_to="2021")
        self.assertEqual(len(articles), 1)                 # the record without a PMID is skipped
        article = articles[0]
        self.assertEqual(article.title, "Lactase persistence in Europe")
        self.assertEqual(article.authors, ["Ann Smith", "Test Consortium"])
        self.assertEqual(article.abstract, "BACKGROUND: Why.\nHow.")   # sections on separate lines
        self.assertEqual((article.doi, article.pmcid, article.journal), ("10.1000/abc", "PMC9", "Journal of Tests"))
        self.assertEqual(article.publication_date, "2020 Mar 04")
        query = parse_qs(urlparse(urls[0]).query)
        self.assertEqual((query["mindate"], query["maxdate"], query["datetype"]), (["2019/01"], ["2021"], ["pdat"]))

    def test_no_hits_bad_xml_and_bad_arguments(self):
        with mock.patch.object(providers, "_request_json", return_value={"esearchresult": {"idlist": []}}):
            self.assertEqual(L.search_pubmed("nothing"), [])
        with mock.patch.object(providers, "_request_json", return_value={"esearchresult": {"idlist": ["1"]}}), \
                mock.patch.object(providers, "_request_text", return_value="<not xml"):
            with self.assertRaises(L.LiteratureSourceError):
                L.search_pubmed("broken")
        with self.assertRaises(ValueError):
            L.search_pubmed("  ")
        with self.assertRaises(ValueError):
            L.search_pubmed("x", max_results=0)


def _article(**fields):
    return L.LiteratureArticle(**fields)


class CombinedSearchTests(unittest.TestCase):
    def test_partial_failure_deduplication_filters_and_canonical_links(self):
        shared = dict(title="Lactase persistence in Europe", doi="10.1000/ABC", authors=["Ann Smith"],
                      journal="J Test", publication_date="2020-03-04", abstract="Text.",
                      keywords=["lactase"], publication_types=["Journal Article"])
        pubmed = [_article(pmid="111", is_open_access=True, **shared)]
        europe = [_article(pmcid="9", is_open_access=True, **shared),
                  _article(title="Old closed paper", doi="10.1000/old", publication_date="1999-01-01",
                           is_open_access=False)]
        query = L.LiteratureQuery("lactase", sources=["pubmed", "europe_pmc", "openalex"],
                                  max_results=5, date_from="2010", open_access_only=True)
        with mock.patch.object(search, "search_pubmed", return_value=pubmed), \
                mock.patch.object(search, "search_europe_pmc", return_value=europe), \
                mock.patch.object(search, "search_openalex", side_effect=L.LiteratureSourceError("HTTP 503")), \
                mock.patch.object(search, "_enrich_missing_metadata", return_value=(0, [])):
            result = L.run_literature_search(query)
        run = result.run
        assert run is not None
        self.assertEqual(run.source_status["openalex"], "FAILED")
        self.assertEqual(run.errors, ["openalex: HTTP 503"])
        self.assertEqual(run.source_status["pubmed"], "OK (1 records)")
        self.assertEqual(run.source_status["europe_pmc"], "OK (1 records)")   # closed, old record filtered
        self.assertEqual(run.total_results_before_deduplication, 2)
        self.assertEqual(len(run.articles), 1)                                 # same DOI merged
        urls = run.articles[0].source_urls
        self.assertIn("https://pubmed.ncbi.nlm.nih.gov/111/", urls)
        self.assertIn("https://pmc.ncbi.nlm.nih.gov/articles/PMC9/", urls)
        self.assertIn("https://doi.org/10.1000/abc", urls)


class EnrichmentTests(unittest.TestCase):
    def test_europe_pmc_pubmed_and_crossref_fill_missing_fields(self):
        by_pmid = _article(title="Record known by PMID", pmid="111")
        by_doi = _article(title="Record known by DOI", doi="10.1000/xyz")

        def fake_json(url, *a, **k):
            if "europepmc" in url and "EXT_ID" in url:
                return {"resultList": {"result": [{"journalTitle": "J Test", "firstPublicationDate": "2020-03-04",
                                                   "abstractText": "From Europe PMC.", "keywordList": {"keyword": ["k1"]},
                                                   "authorString": "Smith A, Jones B", "pmcid": "PMC9",
                                                   "isOpenAccess": "Y"}]}}
            if "europepmc" in url:
                return {"resultList": {"result": []}}
            if "crossref" in url:
                return {"message": {"container-title": ["Crossref Journal"],
                                    "published-online": {"date-parts": [[2019, 7, 1]]},
                                    "abstract": "<jats:p>From Crossref.</jats:p>",
                                    "author": [{"given": "Cy", "family": "Ross"}]}}
            raise AssertionError(url)

        pubmed_xml = PUBMED_XML.replace("Journal of Tests", "PubMed Journal")
        with mock.patch.object(enrichment, "_request_json", side_effect=fake_json), \
                mock.patch.object(enrichment, "_request_text", return_value=pubmed_xml):
            enriched, warnings = L._enrich_missing_metadata([by_pmid, by_doi])
        self.assertEqual(enriched, 2)
        self.assertEqual(warnings, [])
        self.assertEqual((by_pmid.journal, by_pmid.abstract, by_pmid.pmcid, by_pmid.is_open_access),
                         ("J Test", "From Europe PMC.", "PMC9", True))
        self.assertEqual(by_pmid.authors, ["Smith A", "Jones B"])
        self.assertEqual(by_pmid.publication_types, ["Journal Article"])         # from PubMed
        self.assertEqual((by_doi.journal, by_doi.publication_date, by_doi.authors),
                         ("Crossref Journal", "2019-7-1", ["Cy Ross"]))
        self.assertIn("From Crossref.", by_doi.abstract or "")

    def test_failures_become_warnings_and_keep_the_record(self):
        record = _article(title="Record whose lookups fail", pmid="5", doi="10.1000/fail")
        with mock.patch.object(enrichment, "_request_json", side_effect=L.LiteratureSourceError("HTTP 500")), \
                mock.patch.object(enrichment, "_request_text", side_effect=L.LiteratureSourceError("HTTP 500")):
            enriched, warnings = L._enrich_missing_metadata([record])
        self.assertEqual(enriched, 0)
        self.assertEqual(len(warnings), 3)
        self.assertTrue(all("HTTP 500" in warning for warning in warnings))
        self.assertEqual(record.title, "Record whose lookups fail")


class ExportTests(unittest.TestCase):
    def _run(self):
        run = L.LiteratureSearchRun(query=L.LiteratureQuery("lactase"), run_id="LIT-TEST")
        run.articles = [_article(title="Lactase {persistence} & diet", doi="10.1000/abc", pmid="111",
                                 authors=["Ann Smith", "Bo Li"], journal="J Test",
                                 publication_date="2020-03-04", abstract="Text.", keywords=["lactase"])]
        run.source_status = {"pubmed": "OK (1 records)"}
        return run

    def test_all_formats(self):
        run = self._run()
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            paths = {suffix: L.export_results(run, folder / f"out.{suffix}") for suffix in ("json", "csv", "txt", "ris", "bib")}
            data = json.loads(paths["json"].read_text(encoding="utf-8"))
            rows = list(csv.DictReader(paths["csv"].open(encoding="utf-8")))
            text = paths["txt"].read_text(encoding="utf-8")
            ris = paths["ris"].read_text(encoding="utf-8")
            bib = paths["bib"].read_text(encoding="utf-8")
        self.assertEqual(data["run"]["run_id"], "LIT-TEST")
        self.assertEqual(rows[0]["doi"], "10.1000/abc")
        self.assertIn("Lactase", text)
        self.assertIn("TY  - JOUR", ris)
        self.assertIn("AU  - Ann Smith", ris)
        self.assertIn("ER  -", ris)
        self.assertTrue(bib.startswith("@article{"))
        # An unescaped & stops a LaTeX run (fixed in Phase 7); the DOI stays as it is.
        self.assertIn(r"title = {Lactase \{persistence\} \& diet}", bib)
        self.assertIn("doi = {10.1000/abc}", bib)

    def test_bibtex_escapes_every_latex_special_but_keeps_the_doi(self):
        run = self._run()
        run.articles[0].title = "50% of #1 $cost_ratio ~ x^2"
        run.articles[0].doi = "10.1000/a_b%c"
        with tempfile.TemporaryDirectory() as folder:
            bib = L.export_results(run, Path(folder) / "out.bib").read_text(encoding="utf-8")
        self.assertIn(r"50\% of \#1 \$cost\_ratio \textasciitilde{} x\textasciicircum{}2", bib)
        self.assertIn("doi = {10.1000/a_b%c}", bib)

    def test_unsupported_format(self):
        with self.assertRaisesRegex(ValueError, "Supported export formats"):
            L.export_results(self._run(), "out.docx")


if __name__ == "__main__":
    unittest.main()
