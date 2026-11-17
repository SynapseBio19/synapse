"""Terminal dialogs of analyses/literature.py (moved here in P2-2).

The scientific functions stay in ``analyses.literature``; this module draws the
screens, reads the answers and calls them.
"""

from __future__ import annotations

from typing import Any
from typing import Optional
import sys
from analyses.literature import (  # noqa: F401
    LiteratureError,
    LiteratureQuery,
    LiteratureSearchRun,
    MODULE_NAME,
    MODULE_VERSION,
    _normalize_keywords,
    export_results,
    run_literature_search,
    self_test,
)
from ui import style
from ui.style import minor_rule


def print_results(run: LiteratureSearchRun, limit: int = 10) -> None:

    style.header('Literature search results')

    style.field('Query', f"{run.query.query}")

    style.field('Sources', f"{', '.join(run.query.sources)}")

    style.field('Retrieved before deduplication', f"{run.total_results_before_deduplication}")

    style.field('Unique articles', f"{len(run.articles)}")

    print(minor_rule())

    for index, article in enumerate(run.articles[:limit], 1):

        date_flag = " [REVIEW: FUTURE PUBLICATION DATE]" if article.metadata.get("date_validation") == "FUTURE_DATE_REVIEW" else ""
        print(f"{index}. {article.title}{date_flag}")

        print(f"   Relevance: {article.relevance_score if article.relevance_score is not None else 'N/A'} / 100 (text match; not study quality)")
        if article.metadata.get("metadata_completeness") is not None:
            print(f"   Metadata completeness: {article.metadata['metadata_completeness']:.0%} (bibliographic fields only)")
        missing_core = article.metadata.get("missing_core_metadata", [])
        if missing_core:
            print(f"   Missing metadata: {', '.join(missing_core)}")

        if article.authors:
            shown_authors = ", ".join(article.authors[:5])
            if len(article.authors) > 5:
                shown_authors += f" … (showing 5 of {len(article.authors)}; open article detail for the complete list)"
            print(f"   Authors: {shown_authors}")
        else:
            print("   Authors: Not available")

        date_display = article.publication_date or "Not available"
        if article.metadata.get("date_validation") == "FUTURE_DATE_REVIEW":
            date_display += " [REVIEW: FUTURE DATE]"
        print(f"   Journal: {article.journal or 'Not available'} | Date: {date_display}")

        print(f"   DOI: {article.doi or 'Not available'} | PMID: {article.pmid or 'Not available'}")

        print(f"   Abstract: {(article.abstract[:450] + '…') if article.abstract and len(article.abstract) > 450 else article.abstract or 'Not available'}")

        print(f"   URL: {article.source_urls[0] if article.source_urls else 'Not available'}")

        if article.publication_types:
            print(f"   Publication type: {', '.join(article.publication_types[:3])}")
        print(f"   Sources: {', '.join(article.source_ids.keys()) or 'Unknown'}")

        print()

    if run.warnings:

        style.section("WARNINGS", blank_before=False)

        for warning in run.warnings: print(f"- {warning}")

    if run.errors:

        style.section("SOURCE ERRORS", blank_before=bool(run.warnings))

        for error in run.errors: print(f"- {error}")

    print("\nInterpretation: relevance ranking reflects text/query overlap, not scientific quality or truth.")


def print_article_detail(run: LiteratureSearchRun, article_index: int) -> None:

    """Display complete normalized metadata for one result (1-based index)."""

    if article_index < 1 or article_index > len(run.articles):

        style.error("Article number is outside the result range.")

        return

    article = run.articles[article_index - 1]

    style.header('LITERATURE ARTICLE DETAIL')

    style.field('Title', f"{article.title}")

    style.field('Authors', f"{', '.join(article.authors) or 'Not available'}")

    style.field('Journal', f"{article.journal or 'Not available'}")

    date_display = article.publication_date or "Not available"
    if article.metadata.get("date_validation") == "FUTURE_DATE_REVIEW":
        date_display += " [REVIEW: FUTURE DATE]"
    style.field('Publication date', f"{date_display}")
    if article.metadata.get("date_validation_note"):
        print(f"DATE VALIDATION WARNING: {article.metadata['date_validation_note']}")

    style.field('DOI', f"{article.doi or 'Not available'}")

    style.field('PMID', f"{article.pmid or 'Not available'} | PMCID: {article.pmcid or 'Not available'}")

    style.field('Open access', f"{'Yes' if article.is_open_access is True else 'No' if article.is_open_access is False else 'Unknown'}")

    style.field('Citation count', f"{article.citation_count if article.citation_count is not None else 'Not available'}")

    style.field('Relevance', f"{article.relevance_score if article.relevance_score is not None else 'Not scored'} / 100 (text match; not study quality)")
    if article.metadata.get("metadata_completeness") is not None:
        style.field('Metadata completeness', f"{article.metadata['metadata_completeness']:.0%} (bibliographic fields only)")
    if article.metadata.get("missing_core_metadata"):
        style.field('Missing core metadata', f"{', '.join(article.metadata['missing_core_metadata'])}")

    print("Matched terms/components:")

    for key, value in article.relevance_components.items():

        print(f"  {key}: {value}")

    style.field('Publication type(s)', f"{', '.join(article.publication_types) or 'Not available from retrieved metadata'}")
    style.field('Sources', f"{', '.join(article.source_ids) or 'Not available'}")

    print("Links:")

    for url in article.source_urls:

        print(f"  {url}")

    style.section("ABSTRACT")

    print(article.abstract or "Abstract not available from the retrieved sources.")

    style.section("KEYWORDS / SUBJECT TERMS")

    print(", ".join(_normalize_keywords(article.keywords)) or "Not available")

    print()
    style.note("Scientific limitation: metadata and relevance are not an appraisal of study quality or validity.")


def _interactive_menu(context: Optional[dict[str, Any]] = None) -> None:
    context = dict(context or {})
    suggested = str(context.get("suggested_query") or "").strip()

    while True:

        style.header("Literature search")
        style.field("Module", f"{MODULE_NAME} v{MODULE_VERSION}")
        style.field("Sources", "PubMed, Europe PMC, OpenAlex, Crossref")
        if suggested:
            style.field("Suggested query", suggested)
        style.menu([
            ("Search", [("1", "New Literature Search")]),
            ("Configuration", [
                ("2", "Search Configuration (current run parameters are entered per search)"),
            ]),
            ("System", [("3", "Run Self-Test"), ("4", "Data Sources & Scientific Limitations")]),
        ])

        choice = style.ask_choice()

        if choice == "0":

            return

        if choice == "3":

            self_test()

            continue

        if choice == "4":

            style.section("Data sources and scientific limitations")
            style.field("Sources", "PubMed, Europe PMC, OpenAlex, Crossref")
            style.notes([
                "Relevance scores measure query-text matching only; they are not evidence quality.",
                "Source coverage, indexing delay, API availability, and metadata completeness vary.",
                "Abstracts and full text may be unavailable or restricted by copyright/licensing.",
            ])
            style.pause()

            continue

        if choice == "2":

            print("Search options are set when starting each run. Maximum results per source: 200.")

            style.note("Optional environment variables: SYNAPSE_NCBI_EMAIL, SYNAPSE_OPENALEX_EMAIL, SYNAPSE_CROSSREF_EMAIL.")

            continue

        if choice != "1":

            style.invalid_selection()

            continue

        print("Tip: use specific biological terms (gene/protein, organism, disease, pathway) instead of broad single words.")

        print("Example: BRCA1 DNA repair breast cancer Homo sapiens")

        if suggested:
            print(f"Suggestion from the active sequence header: {suggested}")
            search_text = style.ask('Search query [Enter = use suggestion]') or suggested
        else:
            search_text = style.ask('Search query (Boolean terms supported by each source)')

        if not search_text:

            style.error("Query cannot be empty.")

            continue

        max_raw = style.ask('Maximum results per source [30]')

        try:

            max_results = int(max_raw) if max_raw else 30

            source_raw = style.ask('Sources [1=PubMed, 2=Europe PMC, 3=OpenAlex, 4=Crossref, Enter=1,2,3]')

            source_map = {"1": "pubmed", "2": "europe_pmc", "3": "openalex", "4": "crossref"}

            if source_raw:

                selected = [source_map[x.strip()] for x in source_raw.split(",") if x.strip() in source_map]

            else:

                selected = ["pubmed", "europe_pmc", "openalex"]

            if not selected:

                style.error("Select at least one valid source.")

                continue

            date_from = style.ask('From date YYYY or YYYY-MM-DD [none]') or None

            date_to = style.ask('To date YYYY or YYYY-MM-DD [none]') or None

            oa = style.confirm("Open-access only?")

            sort = style.ask('Sort [relevance/date/citations] (relevance)').lower() or "relevance"

            query = LiteratureQuery(search_text, selected, max_results, date_from, date_to, oa, sort,
                                    biological_context={k: v for k, v in context.items() if k == "sequence_header"})

            query.validate()

        except (ValueError, KeyError) as exc:

            style.error(f"Invalid configuration: {exc}")

            continue

        print("\nSearching selected literature sources. A source may fail independently...\n")

        result = run_literature_search(query)

        if not result.run:

            style.error(result.message or "Search did not produce a run.")

            continue

        run = result.run

        if sort == "date":

            run.articles.sort(key=lambda a: a.publication_date or "", reverse=True)

        elif sort == "citations":

            run.articles.sort(key=lambda a: a.citation_count or 0, reverse=True)

        print_results(run)

        while run.articles:

            detail_choice = style.ask('View full article detail? Enter result number, or press Enter to continue')

            if not detail_choice:

                break

            try:

                print_article_detail(run, int(detail_choice))

            except ValueError:

                style.error("Enter a valid result number.")

        export = style.ask('Export results? Enter path with .json/.csv/.txt/.ris/.bib, or press Enter to skip')

        if export:

            try:

                saved = export_results(run, export)

                style.success(f"Exported: {saved.resolve()}")

            except (OSError, ValueError) as exc:

                style.error(f"Export failed: {exc}")


def main(context: Optional[dict[str, Any]] = None) -> int:
    if context is not None:
        _interactive_menu(context)
        return 0

    if len(sys.argv) > 1 and sys.argv[1] in {"--self-test", "--test"}:

        return 0 if self_test() else 1

    if len(sys.argv) > 1 and sys.argv[1] == "--search":

        if len(sys.argv) < 3:

            print('Usage: python literature.py --search "your query" [max_results]')

            return 2

        try:

            max_results = int(sys.argv[3]) if len(sys.argv) > 3 else 20

            query = LiteratureQuery(query=sys.argv[2], max_results=max_results)

            result = run_literature_search(query)

            if result.run:

                print_results(result.run)

                return 0 if result.success else 1

        except (ValueError, LiteratureError) as exc:

            style.error(f"Search failed: {exc}")

            return 1

    _interactive_menu()

    return 0
