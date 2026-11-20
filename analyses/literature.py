"""Synapse Literature Intelligence.

A reproducible, multi-source literature retrieval and analysis module for

scientific workflows. Uses official public APIs and the Python standard library.

Sources: PubMed/NCBI E-utilities, Europe PMC, OpenAlex, Crossref.

Scores indicate query relevance only; they are not evidence quality, truth,

clinical validity, or experimental confidence scores.

"""

# Phase 10: the implementation moved to the internal package
# ``analyses._literature`` (one module per responsibility, see its
# __init__). This module remains the public interface and re-exports
# every name it defined before.

from __future__ import annotations

from analyses._literature.models import (
    DEFAULT_TIMEOUT,
    LiteratureArticle,
    LiteratureError,
    LiteratureQuery,
    LiteratureSearchResult,
    LiteratureSearchRun,
    LiteratureSourceError,
    MAX_RESULTS_LIMIT,
    MODULE_NAME,
    MODULE_VERSION,
    RelevanceComponents,
    USER_AGENT,
)
from analyses._literature.normalize import (
    _apply_date_validation,
    _clean_text,
    _date_bound,
    _date_in_range,
    _deduplicate_author_names,
    _future_date_note,
    _keyword_identity,
    _normalize_keywords,
    _normalize_person_name,
)
from analyses._literature.http import (
    _http_error_is_not_found,
    _LAST_REQUEST_AT,
    MAX_REQUEST_ATTEMPTS,
    MIN_REQUEST_INTERVAL_SECONDS,
    _ncbi_api_key,
    _request_json,
    _request_text,
    _retry_delay,
    RETRY_STATUS_CODES,
    _throttle,
    _url,
)
from analyses._literature.providers import (
    search_crossref,
    search_europe_pmc,
    search_openalex,
    search_pubmed,
)
from analyses._literature.ranking import (
    _coverage,
    deduplicate_articles,
    _is_placeholder_journal,
    _query_phrase,
    rank_articles,
    _tokenize,
)
from analyses._literature.enrichment import (
    _enrich_missing_metadata,
    _journal_from_pubmed_xml,
    _publication_types_from_pubmed_xml,
)
from analyses._literature.search import run_literature_search
from analyses._literature.export import article_to_dict, export_results
from analyses._literature.context import build_sequence_literature_context
from analyses._literature.selftest import self_test


if __name__ == "__main__":

    # The command-line entry point is part of the interface (P2-2).
    from ui.literature_menu import main
    raise SystemExit(main())
