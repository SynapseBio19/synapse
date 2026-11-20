"""The built-in self-test of the literature module."""

from __future__ import annotations

from core.textformat import field_line
from datetime import datetime, timezone
from typing import Any

from analyses._literature.models import LiteratureArticle, LiteratureQuery, MAX_RESULTS_LIMIT
from analyses._literature.normalize import (
    _apply_date_validation,
    _clean_text,
    _date_in_range,
    _deduplicate_author_names,
    _future_date_note,
    _normalize_keywords,
)
from analyses._literature.http import _http_error_is_not_found
from analyses._literature.ranking import (
    deduplicate_articles,
    _is_placeholder_journal,
    rank_articles,
)
from analyses._literature.enrichment import (
    _enrich_missing_metadata,
    _journal_from_pubmed_xml,
    _publication_types_from_pubmed_xml,
)


# The terminal dialogs of this module are in ui/literature_menu.py (P2-2).


def self_test() -> bool:

    """Deterministic offline tests; no network requests are made."""

    tests: list[tuple[str, Any]] = []

    def check(name: str, condition: bool) -> None:

        tests.append((name, condition))

    article = LiteratureArticle(title="CRISPR gene editing", doi="https://doi.org/10.1000/ABC", pmid="123")

    check("doi_normalization", article.normalized_doi() == "10.1000/abc")

    check("identity_keys", "pmid:123" in article.identity_keys())

    merged = deduplicate_articles([article, LiteratureArticle(title="CRISPR gene editing", doi="10.1000/abc", abstract="Example abstract", authors=["A Author"])])

    check("deduplication_by_doi", len(merged) == 1)

    check("deduplication_preserves_abstract", merged[0].abstract == "Example abstract")

    ranked = rank_articles([LiteratureArticle(title="CRISPR gene editing", abstract="CRISPR gene editing methods"), LiteratureArticle(title="Unrelated topic")], "CRISPR gene editing")

    # rank_articles sets every score; None would mean it did not run.
    check("relevance_ranking", (ranked[0].relevance_score or 0.0) >= (ranked[1].relevance_score or 0.0))

    check("exact_phrase_ranking", rank_articles([LiteratureArticle(title="CRISPR gene editing"), LiteratureArticle(title="Editing genes with CRISPR")], "CRISPR gene editing")[0].title == "CRISPR gene editing")

    check("relevance_components_explain_score", "title_coverage" in ranked[0].relevance_components and "matched_query_terms" in ranked[0].relevance_components)

    check("future_date_warning_logic_supported", bool(_future_date_note("2027-01-01", datetime(2026, 10, 26, tzinfo=timezone.utc))))
    check(
        "future_date_parser_supports_pubmed_month_names",
        bool(_future_date_note("2026 Nov 01", datetime(2026, 10, 26, tzinfo=timezone.utc))),
    )
    check(
        "future_date_parser_supports_month_year",
        bool(_future_date_note("November 2026", datetime(2026, 10, 26, tzinfo=timezone.utc))),
    )
    future_test_article = LiteratureArticle(title="Future paper", publication_date="2027-01-01")
    check(
        "future_date_validation_flags_article_for_review",
        _apply_date_validation(future_test_article, today=datetime(2026, 10, 26, tzinfo=timezone.utc))
        and future_test_article.metadata.get("date_validation") == "FUTURE_DATE_REVIEW",
    )
    check(
        "author_deduplication_preserves_full_names_and_initials",
        _deduplicate_author_names(["A. Krull", "Abigail G Krull", "Seth Entriken", "Seth M Entriken"])
        == ["Abigail G Krull", "Seth Entriken", "Seth M Entriken"],
    )
    check("html_abstract_cleanup", _clean_text("<h4>Motivation</h4><p>A &amp; B</p>") == "Motivation\nA & B")
    metadata_merge = deduplicate_articles([
        LiteratureArticle(title="Metadata merge example", doi="10.1000/meta", journal=None, abstract=None),
        LiteratureArticle(title="Metadata merge example", doi="10.1000/meta", journal="Example Journal", abstract="<p>Useful abstract</p>", keywords=["gene"])
    ])
    check("metadata_merge_fills_journal", len(metadata_merge) == 1 and metadata_merge[0].journal == "Example Journal")
    check("metadata_merge_cleans_abstract", len(metadata_merge) == 1 and metadata_merge[0].abstract == "Useful abstract")

    try:

        LiteratureQuery(query="gene", date_from="2026-02-30").validate()

        invalid_calendar_date = False

    except ValueError:

        invalid_calendar_date = True

    check("calendar_date_validation", invalid_calendar_date)

    try:

        LiteratureQuery(query="gene", date_from="2027", date_to="2026").validate()

        reversed_dates = False

    except ValueError:

        reversed_dates = True

    check("date_range_order_validation", reversed_dates)

    try:
        LiteratureQuery(query="gene", date_from="2026-10-15", date_to="2026-10").validate()
        partial_date_range_valid = True
    except ValueError:
        partial_date_range_valid = False
    check("partial_month_end_date_range_validation", partial_date_range_valid)

    try:
        LiteratureQuery(query="gene", date_from="2026-10-20", date_to="2026-10-15").validate()
        reversed_day_range_rejected = False
    except ValueError:
        reversed_day_range_rejected = True
    check("reversed_day_range_validation", reversed_day_range_rejected)

    context_ranked = rank_articles(

        [

            LiteratureArticle(title="CRISPR gene editing in zebrafish"),

            LiteratureArticle(title="CRISPR gene editing in unrelated models"),

        ],

        "CRISPR gene editing",

        {"organism": "zebrafish"},

    )

    check("biological_context_ranking_component", "biological_context_coverage" in context_ranked[0].relevance_components)

    check("biological_context_match_recorded", "zebrafish" in context_ranked[0].relevance_components.get("matched_context_terms", []))

    try:

        LiteratureQuery(query=" ").validate()

        valid_empty_query = False

    except ValueError:

        valid_empty_query = True

    check("query_validation", valid_empty_query)

    try:

        LiteratureQuery(query="gene", max_results=MAX_RESULTS_LIMIT + 1).validate()

        valid_limit = False

    except ValueError:

        valid_limit = True

    check("result_limit_validation", valid_limit)

    empty_enrichment_count, empty_enrichment_warnings = _enrich_missing_metadata([])
    check("metadata_enrichment_empty_input_is_safe", empty_enrichment_count == 0 and not empty_enrichment_warnings)

    check(
        "author_deduplication_initial_variants",
        _deduplicate_author_names(["A. Krull", "Abigail G Krull", "Seth M Entriken", "Seth Entriken"])
        == ["Abigail G Krull", "Seth M Entriken", "Seth Entriken"],
    )
    check(
        "author_deduplication_keeps_distinct_people",
        len(_deduplicate_author_names(["John Smith", "Jane Smith"])) == 2,
    )
    check(
        "author_deduplication_keeps_ambiguous_middle_name_variants",
        _deduplicate_author_names(["Seth Entriken", "Seth M Entriken"])
        == ["Seth Entriken", "Seth M Entriken"],
    )
    check(
        "author_deduplication_expands_initials_only_forms",
        _deduplicate_author_names(["A. G. Krull", "Abigail Grace Krull"])
        == ["Abigail Grace Krull"],
    )
    check(
        "keyword_deduplication_case_insensitive",
        _normalize_keywords(["CRISPR", "crispr", "Zebrafish"]) == ["CRISPR", "Zebrafish"],
    )
    check(
        "keyword_deduplication_unicode_dash_variants",
        _normalize_keywords(["CRISPR-Cas Systems", "CRISPR–Cas Systems", "crIsPr-cas systems"]) == ["CRISPR-Cas Systems"],
    )
    check(
        "keyword_deduplication_preserves_distinct_terms",
        _normalize_keywords(["CRISPR", "CRISPR screening", "Cas9"]) == ["CRISPR", "CRISPR screening", "Cas9"],
    )
    check(
        "keyword_cleanup_handles_empty_and_markup",
        _normalize_keywords(["", "<i>Zebrafish</i>", "zebrafish"]) == ["Zebrafish"],
    )
    check(
        "http_404_classification",
        _http_error_is_not_found(Exception("HTTP 404 from service")),
    )
    check(
        "http_429_not_classified_as_404",
        not _http_error_is_not_found(Exception("HTTP 429 Too Many Requests")),
    )

    same_title_distinct_ids = deduplicate_articles([
        LiteratureArticle(title="A common title", doi="10.1000/one", abstract="First paper"),
        LiteratureArticle(title="A common title", doi="10.1000/two", abstract="Second paper"),
    ])
    check("distinct_doi_records_with_same_title_not_merged", len(same_title_distinct_ids) == 2)

    title_weighted = rank_articles([
        LiteratureArticle(title="CRISPR zebrafish gene editing", abstract="Relevant methods"),
        LiteratureArticle(title="Unrelated title", abstract="CRISPR zebrafish gene editing and detailed methods"),
    ], "CRISPR gene editing zebrafish")
    check("title_weighted_relevance_is_explainable", title_weighted[0].relevance_components.get("ranking_method") is not None)
    check("query_term_coverage_reported", "query_term_coverage" in title_weighted[0].relevance_components)
    check("metadata_completeness_reported", "metadata_completeness" in title_weighted[0].metadata)
    check("platform_name_not_treated_as_journal", _is_placeholder_journal("PubMed"))
    placeholder_journal = rank_articles([LiteratureArticle(title="CRISPR gene editing", authors=["A Author"], journal="PubMed", publication_date="2024", abstract="CRISPR gene editing", doi="10.1000/test", keywords=["CRISPR"])], "CRISPR gene editing")[0]
    check("placeholder_journal_does_not_inflate_completeness", placeholder_journal.journal is None and placeholder_journal.metadata["metadata_completeness"] < 1.0)
    journal_merge = deduplicate_articles([
        LiteratureArticle(title="Journal merge", doi="10.1000/journal", journal="PubMed"),
        LiteratureArticle(title="Journal merge", doi="10.1000/journal", journal="Zebrafish")
    ])
    check("real_journal_replaces_platform_label", len(journal_merge) == 1 and journal_merge[0].journal == "Zebrafish")
    check("citation_count_does_not_override_relevance", rank_articles([
        LiteratureArticle(title="Unrelated title", citation_count=99999),
        LiteratureArticle(title="CRISPR gene editing zebrafish", citation_count=0),
    ], "CRISPR gene editing zebrafish")[0].title == "CRISPR gene editing zebrafish")

    sample_pubmed_xml = """<PubmedArticleSet><PubmedArticle><MedlineCitation><Article><Journal><Title>Zebrafish</Title></Journal></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"""
    check("pubmed_xml_journal_extraction", _journal_from_pubmed_xml(sample_pubmed_xml) == "Zebrafish")
    check("pubmed_xml_invalid_is_safe", _journal_from_pubmed_xml("<broken") is None)
    check("pubmed_xml_placeholder_journal_rejected", _journal_from_pubmed_xml(
        "<PubmedArticleSet><PubmedArticle><MedlineCitation><Article><Journal><Title>PubMed</Title></Journal></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"
    ) is None)
    pubtype_xml = """<PubmedArticleSet><PubmedArticle><MedlineCitation><Article><PublicationTypeList><PublicationType UI="D016428">Journal Article</PublicationType><PublicationType UI="D016454">Review</PublicationType><PublicationType UI="D016428">Journal Article</PublicationType></PublicationTypeList></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"""
    check("pubmed_publication_type_extraction_and_deduplication",
          _publication_types_from_pubmed_xml(pubtype_xml) == ["Journal Article", "Review"])
    check("pubmed_publication_type_parser_handles_invalid_xml",
          _publication_types_from_pubmed_xml("<broken") == [])
    link_article = LiteratureArticle(title="Canonical links", pmid="12345", pmcid="PMC67890", doi="10.1000/test")
    canonical_urls = [f"https://pubmed.ncbi.nlm.nih.gov/{link_article.pmid}/",
                      f"https://pmc.ncbi.nlm.nih.gov/articles/{link_article.pmcid}/",
                      f"https://doi.org/{link_article.normalized_doi()}"]
    check("canonical_source_links_are_constructible",
          canonical_urls == ["https://pubmed.ncbi.nlm.nih.gov/12345/",
                             "https://pmc.ncbi.nlm.nih.gov/articles/PMC67890/",
                             "https://doi.org/10.1000/test"])
    normalized_article = LiteratureArticle(
        title="Keyword finalization", keywords=["CRISPR-Cas Systems", "CRISPR–Cas Systems", "Zebrafish"],
        publication_types=["Journal Article", "Journal Article", " Review "],
    )
    normalized_article.keywords = _normalize_keywords(normalized_article.keywords)
    normalized_article.publication_types = list(dict.fromkeys(_clean_text(x) for x in normalized_article.publication_types if _clean_text(x)))
    check("final_keyword_normalization_removes_unicode_duplicates", normalized_article.keywords == ["CRISPR-Cas Systems", "Zebrafish"])
    check("publication_types_are_cleaned_and_deduplicated", normalized_article.publication_types == ["Journal Article", "Review"])

    # Regression tests for date formats and whitespace normalization that the
    # earlier suite did not exercise.
    date_query = LiteratureQuery(query="zebrafish", date_from="2026-10-01", date_to="2026-10-15")
    check("date_filter_iso_day_inside_range", _date_in_range("2026-10-08", date_query))
    check("date_filter_pubmed_month_name_inside_range", _date_in_range("2026 Oct 08", date_query))
    check("date_filter_month_year_overlaps_range", _date_in_range("October 2026", date_query))
    check("date_filter_excludes_later_pubmed_day", not _date_in_range("2026 Oct 20", date_query))
    check("date_filter_year_only_overlaps_range", _date_in_range("2026", date_query))
    check("date_filter_excludes_prior_year", not _date_in_range("2025", date_query))
    check("date_filter_unparseable_is_not_silently_included", not _date_in_range("date unavailable", date_query))
    check("author_normalization_collapses_repeated_whitespace", _deduplicate_author_names(["  Abigail   Grace   Krull  ", "Abigail Grace Krull"]) == ["Abigail Grace Krull"])
    check("author_normalization_handles_tabs", _deduplicate_author_names(["Abigail\tGrace Krull", "Abigail Grace Krull"]) == ["Abigail Grace Krull"])
    check("future_date_not_flagged_for_past_pubmed_month", _future_date_note("2026 Oct 08", datetime(2026, 10, 27, tzinfo=timezone.utc)) is None)
    check("future_date_flags_future_pubmed_month", bool(_future_date_note("2026 Nov 01", datetime(2026, 10, 27, tzinfo=timezone.utc))))

    invalid_day_query = LiteratureQuery(query="zebrafish", date_from="2026-02-01", date_to="2026-02-28")
    check("date_filter_rejects_impossible_day", not _date_in_range("2026-02-31", invalid_day_query))
    check("date_filter_rejects_impossible_month", not _date_in_range("2026-13", invalid_day_query))
    check("future_date_does_not_flag_impossible_day", _future_date_note("2026-11-31", datetime(2026, 10, 27, tzinfo=timezone.utc)) is None)

    for name, passed in tests:

        print(f"[{'PASS' if passed else 'FAIL'}] {name}")

    print(field_line('Self-test', f"{sum(ok for _, ok in tests)}/{len(tests)} passed"))

    return all(ok for _, ok in tests)
