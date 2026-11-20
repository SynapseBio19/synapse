"""The combined search over the selected providers."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

from analyses._literature.models import (
    LiteratureArticle,
    LiteratureQuery,
    LiteratureSearchResult,
    LiteratureSearchRun,
    LiteratureSourceError,
)
from analyses._literature.normalize import (
    _apply_date_validation,
    _clean_text,
    _date_in_range,
    _deduplicate_author_names,
    _normalize_keywords,
)
from analyses._literature.providers import (
    search_crossref,
    search_europe_pmc,
    search_openalex,
    search_pubmed,
)
from analyses._literature.ranking import (
    deduplicate_articles,
    _is_placeholder_journal,
    rank_articles,
)
from analyses._literature.enrichment import _enrich_missing_metadata


def run_literature_search(query: LiteratureQuery, email: Optional[str] = None) -> LiteratureSearchResult:

    """Execute selected sources independently, preserving partial results if one fails."""

    query.validate()

    run = LiteratureSearchRun(query=query, run_id=datetime.now(timezone.utc).strftime("LIT-%Y%m%dT%H%M%S-%fZ"))

    collectors = {

        "pubmed": lambda: search_pubmed(query.query, query.max_results, query.date_from, query.date_to, email),

        "europe_pmc": lambda: search_europe_pmc(query.query, query.max_results, query.date_from, query.date_to, query.open_access_only),

        "openalex": lambda: search_openalex(query.query, query.max_results, query.date_from, query.date_to, email),

        "crossref": lambda: search_crossref(query.query, query.max_results, query.date_from, query.date_to, email),

    }

    collected: list[LiteratureArticle] = []

    for source in query.sources:

        try:

            records = collectors[source]()

            if query.open_access_only:

                records = [r for r in records if r.is_open_access is True]

            if query.date_from or query.date_to:

                records = [r for r in records if _date_in_range(r.publication_date, query)]

            collected.extend(records)

            run.source_status[source] = f"OK ({len(records)} records)"

        except (LiteratureSourceError, ValueError) as exc:

            run.source_status[source] = "FAILED"

            run.errors.append(f"{source}: {exc}")

    run.total_results_before_deduplication = len(collected)

    run.articles = deduplicate_articles(collected)

    # Try to fill missing bibliographic fields from independent services.
    # This is bounded to avoid excessive API traffic and never invalidates results.
    enriched_count, enrichment_warnings = _enrich_missing_metadata(run.articles)
    if enriched_count:
        run.warnings.append(
            f"Metadata enrichment filled missing fields for {enriched_count} record(s) "
            "using Europe PMC and/or Crossref."
        )
    run.warnings.extend(enrichment_warnings)

    # Final, source-independent normalization pass. This also protects against source
    # adapters or enrichment paths that introduce duplicate labels after initial merging.
    for article in run.articles:
        article.authors = _deduplicate_author_names(article.authors)
        article.keywords = _normalize_keywords(article.keywords)
        article.publication_types = list(dict.fromkeys(
            _clean_text(item) for item in article.publication_types if _clean_text(item)
        ))
        if _is_placeholder_journal(article.journal):
            article.journal = None
        # Ensure users can reach canonical records even when a source adapter
        # only retained a DOI URL or omitted its source URL during merging.
        canonical_urls = []
        if article.pmid:
            canonical_urls.append(f"https://pubmed.ncbi.nlm.nih.gov/{article.pmid.strip()}/")
        if article.pmcid:
            pmcid = article.pmcid.strip()
            if not pmcid.upper().startswith("PMC"):
                pmcid = "PMC" + pmcid
            canonical_urls.append(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/")
        doi = article.normalized_doi()
        if doi:
            canonical_urls.append(f"https://doi.org/{doi}")
        article.source_urls = list(dict.fromkeys(article.source_urls + canonical_urls))

    # Recalculate date flags after metadata enrichment may have supplied dates.
    run.articles = rank_articles(run.articles, query.query, query.biological_context)

    # Flag suspicious future dates without silently changing source-provided values.
    future_dated = [article for article in run.articles if _apply_date_validation(article)]
    if future_dated:
        # Name the exact records so users can investigate the warning without
        # having to guess which result triggered it. Keep source dates unchanged.
        future_details = []
        for article in future_dated:
            title = _clean_text(article.title) or "Untitled record"
            date_value = _clean_text(article.publication_date) or "Unknown date"
            future_details.append(f"{title} [{date_value}]")
        details_text = " | ".join(future_details)
        run.warnings.append(
            f"{len(future_dated)} record(s) report a publication date later than today: "
            f"{details_text}. Each record is flagged [REVIEW: FUTURE DATE]; "
            "source-provided dates were not changed. Check the publisher record and source metadata."
        )
    if not run.articles and run.errors:

        run.warnings.append("No results were retrieved; all selected sources may have failed. Check network access and source status.")

    if any(status == "FAILED" for status in run.source_status.values()):

        run.warnings.append("One or more sources failed. Results may be incomplete; inspect source status and errors.")

    if not any(article.abstract for article in run.articles):

        run.warnings.append("No abstracts were available in the retrieved records. This does not mean the articles lack abstracts at the publisher.")

    run.finish()

    return LiteratureSearchResult(success=bool(run.articles), run=run,

                                  message=f"Retrieved {len(run.articles)} unique articles from {len(query.sources)} selected source(s).")
