"""Split plan for analyses/literature.py (Phase 10)."""

SOURCE = "analyses/literature.py"
PACKAGE = "analyses._literature"
PACKAGE_DIR = "analyses/_literature"
PACKAGE_DOC = """Implementation of :mod:`analyses.literature` (split out in Phase 10).

Import from ``analyses.literature``; it re-exports every name defined here
and is the only public interface. The modules, in dependency order:

- ``models``: module constants, errors, the query, article and run records
- ``normalize``: text, name, keyword and date normalization and checks
- ``http``: the HTTP layer with per-host throttling and bounded retry
- ``providers``: the PubMed, Europe PMC, OpenAlex and Crossref searches
- ``ranking``: de-duplication and the transparent relevance ranking
- ``enrichment``: filling missing metadata from Europe PMC, PubMed, Crossref
- ``search``: the combined search ``run_literature_search``
- ``export``: JSON, CSV, text, RIS and BibTeX exports
- ``context``: literature context of a selected sequence or ORF
- ``selftest``: the built-in self-test

Tests that replace a helper must replace it in the module that calls it
(for example ``analyses._literature.providers._request_json``); a stand-in
on ``analyses.literature`` would not reach the moved code.
"""

MODULES = {
    "models": ("Module constants, errors, and the query, article and run records.", [
        "RelevanceComponents", "MODULE_NAME", "MODULE_VERSION", "USER_AGENT", "DEFAULT_TIMEOUT",
        "MAX_RESULTS_LIMIT", "LiteratureError", "LiteratureSourceError", "LiteratureQuery",
        "LiteratureArticle", "LiteratureSearchRun", "LiteratureSearchResult",
    ]),
    "normalize": ("Text, person-name, keyword and date normalization and checks.", [
        "_date_bound", "_clean_text", "_normalize_person_name", "_deduplicate_author_names",
        "_keyword_identity", "_normalize_keywords", "_date_in_range", "_future_date_note",
        "_apply_date_validation",
    ]),
    "http": ("The HTTP layer: per-host throttling, bounded retry, JSON and URL helpers.", [
        "MIN_REQUEST_INTERVAL_SECONDS", "RETRY_STATUS_CODES", "MAX_REQUEST_ATTEMPTS", "_LAST_REQUEST_AT",
        "_ncbi_api_key", "_throttle", "_retry_delay", "_request_text", "_request_json", "_url",
        "_http_error_is_not_found",
    ]),
    "providers": ("The PubMed, Europe PMC, OpenAlex and Crossref searches.", [
        "search_pubmed", "search_europe_pmc", "search_openalex", "search_crossref",
    ]),
    "ranking": ("De-duplication of records and the transparent relevance ranking.", [
        "_is_placeholder_journal", "deduplicate_articles", "_tokenize", "_coverage", "_query_phrase",
        "rank_articles",
    ]),
    "enrichment": ("Filling missing metadata from Europe PMC, PubMed and Crossref.", [
        "_journal_from_pubmed_xml", "_publication_types_from_pubmed_xml", "_enrich_missing_metadata",
    ]),
    "search": ("The combined search over the selected providers.", [
        "run_literature_search",
    ]),
    "export": ("JSON, CSV, text, RIS and BibTeX exports.", [
        "article_to_dict", "export_results",
    ]),
    "context": ("Literature context of a selected sequence or ORF.", [
        "build_sequence_literature_context",
    ]),
    "selftest": ("The built-in self-test of the literature module.", [
        "self_test",
    ]),
}

FACADE_NAMES = set()
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 10: the implementation moved to the internal package\n"
    "# ``analyses._literature`` (one module per responsibility, see its\n"
    "# __init__). This module remains the public interface and re-exports\n"
    "# every name it defined before.\n"
)
