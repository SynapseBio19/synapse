"""Module constants, errors, and the query, article and run records."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Optional, TypedDict


class RelevanceComponents(TypedDict, total=False):
    """How a relevance score was composed (rank_articles); keys appear as computed."""

    title_coverage: float
    abstract_coverage: float
    keyword_coverage: float
    query_term_coverage: float
    exact_phrase_bonus: float
    biological_context_coverage: float
    matched_query_terms: list[str]
    missing_query_terms: list[str]
    matched_context_terms: list[str]
    ranking_method: str


MODULE_NAME = "Synapse Literature Intelligence"

MODULE_VERSION = "2.1.1"

USER_AGENT = "Synapse-Literature-Intelligence/2.1.1 (scientific literature retrieval)"

DEFAULT_TIMEOUT = 25

MAX_RESULTS_LIMIT = 200

class LiteratureError(RuntimeError):

    """Base exception for literature retrieval and processing errors."""

class LiteratureSourceError(LiteratureError):

    """Raised when a source request or response fails."""

@dataclass

class LiteratureQuery:

    query: str

    sources: list[str] = field(default_factory=lambda: ["pubmed", "europe_pmc", "openalex"])

    max_results: int = 30

    date_from: Optional[str] = None

    date_to: Optional[str] = None

    open_access_only: bool = False

    sort_by: str = "relevance"  # relevance, date, citations

    biological_context: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:

        if not self.query or not self.query.strip():

            raise ValueError("Search query cannot be empty.")

        if not 1 <= int(self.max_results) <= MAX_RESULTS_LIMIT:

            raise ValueError(f"max_results must be between 1 and {MAX_RESULTS_LIMIT}.")

        allowed = {"pubmed", "europe_pmc", "openalex", "crossref"}

        unknown = set(self.sources) - allowed

        if unknown:

            raise ValueError(f"Unsupported source(s): {', '.join(sorted(unknown))}")

        if not self.sources:

            raise ValueError("Select at least one literature source.")

        if self.sort_by not in {"relevance", "date", "citations"}:

            raise ValueError("sort_by must be relevance, date, or citations.")

        parsed_dates: dict[str, datetime] = {}

        for label, value in (("date_from", self.date_from), ("date_to", self.date_to)):

            if not value:

                continue

            if not re.fullmatch(r"\d{4}(?:-\d{2}(?:-\d{2})?)?", value):

                raise ValueError(f"{label} must use YYYY, YYYY-MM, or YYYY-MM-DD format.")

            parts = [int(part) for part in value.split("-")]

            try:

                parsed = datetime(

                    parts[0],

                    parts[1] if len(parts) > 1 else 1,

                    parts[2] if len(parts) > 2 else 1,

                    tzinfo=timezone.utc,

                )

            except ValueError as exc:

                raise ValueError(f"{label} is not a valid calendar date: {value}") from exc

            parsed_dates[label] = parsed

        # Both bounds are in parsed_dates exactly when both strings are set.
        if self.date_from and self.date_to:
            # Interpret partial dates as ranges: YYYY means the whole year,
            # and YYYY-MM means the whole month. This avoids rejecting a valid
            # range such as 2026-10-15 through 2026-10.
            from_parts = [int(part) for part in self.date_from.split("-")]
            to_parts = [int(part) for part in self.date_to.split("-")]

            from_date = datetime(
                from_parts[0],
                from_parts[1] if len(from_parts) > 1 else 1,
                from_parts[2] if len(from_parts) > 2 else 1,
                tzinfo=timezone.utc,
            )

            to_year = to_parts[0]
            if len(to_parts) == 1:
                to_month, to_day = 12, 31
            elif len(to_parts) == 2:
                import calendar
                to_month = to_parts[1]
                to_day = calendar.monthrange(to_year, to_month)[1]
            else:
                to_month, to_day = to_parts[1], to_parts[2]

            to_date = datetime(to_year, to_month, to_day, tzinfo=timezone.utc)

            if from_date > to_date:
                raise ValueError("date_from must not be later than date_to.")

@dataclass

class LiteratureArticle:

    title: str

    article_id: str = ""

    doi: Optional[str] = None

    pmid: Optional[str] = None

    pmcid: Optional[str] = None

    authors: list[str] = field(default_factory=list)

    journal: Optional[str] = None

    publication_date: Optional[str] = None

    abstract: Optional[str] = None

    keywords: list[str] = field(default_factory=list)

    publication_types: list[str] = field(default_factory=list)

    source_ids: dict[str, str] = field(default_factory=dict)

    source_urls: list[str] = field(default_factory=list)

    is_open_access: Optional[bool] = None

    citation_count: Optional[int] = None

    relevance_score: Optional[float] = None

    relevance_components: RelevanceComponents = field(default_factory=lambda: RelevanceComponents())

    metadata: dict[str, Any] = field(default_factory=dict)

    def normalized_doi(self) -> Optional[str]:

        if not self.doi:

            return None

        value = self.doi.strip().lower()

        value = re.sub(r"^(?:https?://)?(?:dx\.)?doi\.org/", "", value, flags=re.I)
        value = re.sub(r"^doi:\s*", "", value, flags=re.I)

        return value.strip().rstrip(".,; ") or None

    def identity_keys(self) -> list[str]:

        keys: list[str] = []

        doi = self.normalized_doi()

        if doi:

            keys.append("doi:" + doi)

        if self.pmid:

            keys.append("pmid:" + self.pmid.strip())

        if self.pmcid:

            keys.append("pmcid:" + self.pmcid.strip().upper())

        if self.title:

            normalized_title = re.sub(r"[^a-z0-9]+", " ", self.title.lower()).strip()

            if normalized_title:

                keys.append("title:" + normalized_title)

        return keys

@dataclass

class LiteratureSearchRun:

    query: LiteratureQuery

    run_id: str = ""

    started_at_utc: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    completed_at_utc: Optional[str] = None

    articles: list[LiteratureArticle] = field(default_factory=list)

    source_status: dict[str, str] = field(default_factory=dict)

    warnings: list[str] = field(default_factory=list)

    errors: list[str] = field(default_factory=list)

    total_results_before_deduplication: int = 0

    total_results_after_deduplication: int = 0

    def finish(self) -> None:

        self.completed_at_utc = datetime.now(timezone.utc).isoformat()

        self.total_results_after_deduplication = len(self.articles)

@dataclass

class LiteratureSearchResult:

    success: bool

    run: Optional[LiteratureSearchRun] = None

    message: str = ""
