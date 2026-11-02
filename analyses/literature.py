"""Synapse Literature Intelligence.

A reproducible, multi-source literature retrieval and analysis module for

scientific workflows. Uses official public APIs and the Python standard library.

Sources: PubMed/NCBI E-utilities, Europe PMC, OpenAlex, Crossref.

Scores indicate query relevance only; they are not evidence quality, truth,

clinical validity, or experimental confidence scores.

"""

from __future__ import annotations

import csv
import html

import json

import os

import re

import sys

import time

import urllib.error

import urllib.parse
import unicodedata

import urllib.request

import xml.etree.ElementTree as ET

from dataclasses import asdict, dataclass, field

from datetime import datetime, timezone

from pathlib import Path

from typing import Any, Optional
from ui.style import major_rule, minor_rule
from ui import style

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

        if "date_from" in parsed_dates and "date_to" in parsed_dates:
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

    relevance_components: dict[str, float] = field(default_factory=dict)

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

# Minimum interval between requests to the same host (seconds). NCBI
# E-utilities allow 3 requests/s without an API key and 10/s with one
# (SYNAPSE_NCBI_API_KEY); the other services ask for moderate, polite use.
MIN_REQUEST_INTERVAL_SECONDS = {
    "eutils.ncbi.nlm.nih.gov": 0.34,
    "www.ebi.ac.uk": 0.2,
    "api.openalex.org": 0.2,
    "api.crossref.org": 0.2,
}
RETRY_STATUS_CODES = {429, 500, 502, 503, 504}
MAX_REQUEST_ATTEMPTS = 3
_LAST_REQUEST_AT: dict[str, float] = {}


def _ncbi_api_key() -> Optional[str]:
    key = os.getenv("SYNAPSE_NCBI_API_KEY", "").strip()
    return key or None


def _throttle(host: str) -> None:
    interval = MIN_REQUEST_INTERVAL_SECONDS.get(host, 0.2)
    if host == "eutils.ncbi.nlm.nih.gov" and _ncbi_api_key():
        interval = 0.1
    last = _LAST_REQUEST_AT.get(host)
    if last is not None:
        wait = interval - (time.monotonic() - last)
        if wait > 0:
            time.sleep(wait)
    _LAST_REQUEST_AT[host] = time.monotonic()


def _retry_delay(exc: urllib.error.HTTPError, attempt: int) -> float:
    retry_after = exc.headers.get("Retry-After") if exc.headers else None
    if retry_after and retry_after.strip().isdigit():
        return min(30.0, float(retry_after))
    return min(30.0, 2.0 ** attempt)


def _request_text(url: str, timeout: int = DEFAULT_TIMEOUT, headers: Optional[dict[str, str]] = None) -> str:
    """GET a URL with per-host throttling and bounded retry on 429/5xx."""

    request_headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}

    if headers:

        request_headers.update(headers)

    host = urllib.parse.urlsplit(url).netloc

    for attempt in range(1, MAX_REQUEST_ATTEMPTS + 1):

        _throttle(host)

        request = urllib.request.Request(url, headers=request_headers)

        try:

            with urllib.request.urlopen(request, timeout=timeout) as response:

                charset = response.headers.get_content_charset() or "utf-8"

                return response.read().decode(charset, errors="replace")

        except urllib.error.HTTPError as exc:

            if exc.code in RETRY_STATUS_CODES and attempt < MAX_REQUEST_ATTEMPTS:

                time.sleep(_retry_delay(exc, attempt))

                continue

            raise LiteratureSourceError(

                f"HTTP {exc.code} from {host} after {attempt} attempt(s); "

                "the source may be rate-limiting requests or temporarily unavailable."

            ) from exc

        except (urllib.error.URLError, TimeoutError, OSError) as exc:

            raise LiteratureSourceError(f"Request failed for {host}: {exc}") from exc

    raise LiteratureSourceError(f"Request failed for {host}.")


def _date_bound(value: Optional[str], upper: bool) -> Optional[str]:
    """Expand YYYY / YYYY-MM / YYYY-MM-DD to a full ISO date (range start or end)."""
    if not value:
        return None
    import calendar
    parts = [int(p) for p in value.split("-")]
    year = parts[0]
    if len(parts) == 1:
        month, day = (12, 31) if upper else (1, 1)
    elif len(parts) == 2:
        month = parts[1]
        day = calendar.monthrange(year, month)[1] if upper else 1
    else:
        month, day = parts[1], parts[2]
    return f"{year:04d}-{month:02d}-{day:02d}"

def _request_json(url: str, timeout: int = DEFAULT_TIMEOUT, headers: Optional[dict[str, str]] = None) -> dict[str, Any]:

    try:

        payload = json.loads(_request_text(url, timeout=timeout, headers=headers))

    except json.JSONDecodeError as exc:

        raise LiteratureSourceError("The source returned invalid JSON.") from exc

    if not isinstance(payload, dict):

        raise LiteratureSourceError("The source returned an unexpected JSON structure.")

    return payload

def _url(base: str, params: dict[str, Any]) -> str:

    return base + "?" + urllib.parse.urlencode({k: v for k, v in params.items() if v is not None})

def _clean_text(value: Any) -> str:
    """Normalize text from APIs, removing HTML markup and decoding entities."""
    if value is None:
        return ""
    if isinstance(value, list):
        return " ".join(_clean_text(item) for item in value).strip()
    text = str(value)
    text = re.sub(r"(?is)<(script|style)\b[^>]*>.*?</\1\s*>", " ", text)
    text = re.sub(r"(?i)<br\s*/?>|</(?:p|div|h[1-6]|li|section|article)\s*>", "\n", text)
    text = re.sub(r"(?s)<[^>]+>", " ", text)
    text = html.unescape(text)
    text = text.replace("\r", "\n")
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r" *\n *", "\n", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()



def _normalize_person_name(name: str) -> str:
    """Normalize a person name for conservative duplicate detection."""
    value = unicodedata.normalize("NFKD", str(name or ""))
    value = "".join(ch for ch in value if not unicodedata.combining(ch))
    value = value.casefold().replace("&", " and ")
    value = re.sub(r"\b(?:jr|sr|ii|iii|iv)\b\.?", " ", value)
    value = re.sub(r"[^\w\s]", " ", value)
    return " ".join(value.split())


def _deduplicate_author_names(names: list[str]) -> list[str]:
    """Deduplicate exact labels and clear initials-only variants conservatively.

    A full given name with an omitted middle name is not treated as equivalent
    to a longer full-name form: that match is ambiguous without an external ID.
    """
    output: list[str] = []
    keys: list[tuple[str, tuple[str, ...], str]] = []

    def parts(name: str) -> tuple[str, tuple[str, ...], str]:
        cleaned = re.sub(r"\s+", " ", str(name or "")).strip(" ,;")
        if not cleaned:
            return "", (), ""
        if "," in cleaned:
            family, given = [part.strip() for part in cleaned.split(",", 1)]
        else:
            tokens = cleaned.split()
            family, given = tokens[-1], " ".join(tokens[:-1])
        family_key = _normalize_person_name(family)
        given_tokens = tuple(_normalize_person_name(given).split())
        normalized = " ".join(given_tokens)
        return family_key, given_tokens, normalized

    def initials_only(tokens: tuple[str, ...]) -> bool:
        return bool(tokens) and all(len(token) == 1 for token in tokens)

    def compatible_initials(left: tuple[str, ...], right: tuple[str, ...]) -> bool:
        # Only allow expansion when one entire given-name representation is
        # initials-only, e.g. "A. G. Krull" vs "Abigail Grace Krull".
        if not (initials_only(left) or initials_only(right)):
            return False
        short, long = (left, right) if initials_only(left) else (right, left)
        if len(short) > len(long):
            return False
        return all(len(token) == 1 and token[0] == long[index][0]
                   for index, token in enumerate(short))

    for raw in names or []:
        name = re.sub(r"\s+", " ", str(raw or "")).strip(" ,;")
        if not name:
            continue
        family_key, given_tokens, normalized_given = parts(name)
        if not family_key:
            continue
        duplicate_index = None
        for index, (existing_family, existing_given_tokens, existing_normalized) in enumerate(keys):
            if family_key != existing_family:
                continue
            if normalized_given == existing_normalized:
                duplicate_index = index
                break
            if compatible_initials(given_tokens, existing_given_tokens):
                duplicate_index = index
                break

        if duplicate_index is None:
            output.append(name)
            keys.append((family_key, given_tokens, normalized_given))
        elif len(name) > len(output[duplicate_index]):
            output[duplicate_index] = name
            keys[duplicate_index] = parts(name)

    return output


def _keyword_identity(value: str) -> str:
    """Create a Unicode- and punctuation-tolerant identity key for subject terms."""
    normalized = unicodedata.normalize("NFKC", _clean_text(value or "")).casefold()
    # Treat typographic hyphens/dashes and punctuation as equivalent separators.
    normalized = normalized.replace("&", " and ")
    normalized = re.sub(r"[‐‑‒–—−-]+", " ", normalized)
    normalized = re.sub(r"[^\w\s]", " ", normalized, flags=re.UNICODE)
    return " ".join(normalized.split())


def _normalize_keywords(keywords: list[str]) -> list[str]:
    """Remove duplicate subject terms across source variants, preserving first display form."""
    output: list[str] = []
    seen: set[str] = set()
    for keyword in keywords or []:
        value = _clean_text(str(keyword or "")).strip(" ;,\t\n")
        key = _keyword_identity(value)
        if value and key and key not in seen:
            seen.add(key)
            output.append(value)
    return output


def _http_error_is_not_found(exc: Exception) -> bool:
    """Return True for a definite HTTP 404 response."""
    return "HTTP 404" in str(exc) or "status code 404" in str(exc).casefold()

def _date_in_range(date_value: Optional[str], query: LiteratureQuery) -> bool:
    """Apply inclusive date bounds to ISO and common bibliographic date formats.

    Publication dates often have reduced precision (year or year/month) or use
    PubMed's ``YYYY Mon DD`` format. A partial publication date is treated as an
    interval, so it is included when that interval overlaps the requested range.
    Unparseable dates are retained only when no date filter was requested.
    """
    if not date_value:
        return not (query.date_from or query.date_to)

    value = _clean_text(date_value).strip()
    year_match = re.search(r"(?<!\d)((?:19|20)\d{2})(?!\d)", value)
    if not year_match:
        return not (query.date_from or query.date_to)

    year = int(year_match.group(1))
    month = day = None
    # Numeric forms: YYYY-MM-DD, YYYY/MM/DD, YYYY MM DD.
    numeric = re.search(r"(?:19|20)\d{2}[- /](\d{1,2})(?:[- /](\d{1,2}))?", value)
    if numeric:
        candidate_month = int(numeric.group(1))
        candidate_day = int(numeric.group(2)) if numeric.group(2) else None
        if not 1 <= candidate_month <= 12:
            return not (query.date_from or query.date_to)
        month = candidate_month
        if candidate_day is not None:
            if not 1 <= candidate_day <= 31:
                return not (query.date_from or query.date_to)
            day = candidate_day
    else:
        month_names = {
            "jan": 1, "january": 1, "feb": 2, "february": 2,
            "mar": 3, "march": 3, "apr": 4, "april": 4,
            "may": 5, "jun": 6, "june": 6, "jul": 7, "july": 7,
            "aug": 8, "august": 8, "sep": 9, "sept": 9, "september": 9,
            "oct": 10, "october": 10, "nov": 11, "november": 11,
            "dec": 12, "december": 12,
        }
        month_match = re.search(r"\b(" + "|".join(month_names) + r")\b", value, re.I)
        if month_match:
            month = month_names[month_match.group(1).casefold()]
            remainder = value.replace(year_match.group(1), " ", 1)
            remainder = re.sub(re.escape(month_match.group(1)), " ", remainder, count=1, flags=re.I)
            day_match = re.search(r"\b(\d{1,2})\b", remainder)
            if day_match and 1 <= int(day_match.group(1)) <= 31:
                day = int(day_match.group(1))

    import calendar
    if month is not None:
        try:
            last_day = calendar.monthrange(year, month)[1]
        except ValueError:
            return not (query.date_from or query.date_to)
        if day is not None:
            if day > last_day:
                # An impossible source date is unparseable, not a valid month-level date.
                return not (query.date_from or query.date_to)
            article_low = article_high = (year, month, day)
        else:
            article_low, article_high = (year, month, 1), (year, month, last_day)
    else:
        article_low, article_high = (year, 1, 1), (year, 12, 31)

    def lower_bound(raw: str) -> tuple[int, int, int]:
        parts = [int(part) for part in raw.split("-")]
        return parts[0], parts[1] if len(parts) > 1 else 1, parts[2] if len(parts) > 2 else 1

    def upper_bound(raw: str) -> tuple[int, int, int]:
        parts = [int(part) for part in raw.split("-")]
        if len(parts) == 1:
            return parts[0], 12, 31
        if len(parts) == 2:
            return parts[0], parts[1], calendar.monthrange(parts[0], parts[1])[1]
        return parts[0], parts[1], parts[2]

    if query.date_from and article_high < lower_bound(query.date_from):
        return False
    if query.date_to and article_low > upper_bound(query.date_to):
        return False
    return True

def _future_date_note(publication_date: Optional[str], today: Optional[datetime] = None) -> Optional[str]:
    """Return a review note when a source-supplied publication date is in the future."""
    if not publication_date:
        return None
    today = today or datetime.now(timezone.utc)
    value = _clean_text(publication_date).strip()
    year_match = re.search(r"(?:19|20)\d{2}", value)
    if not year_match:
        return None
    year = int(year_match.group(0))
    month = day = None

    numeric = re.search(r"(?:19|20)\d{2}[- /](\d{1,2})(?:[- /](\d{1,2}))?", value)
    if numeric:
        month = int(numeric.group(1))
        day = int(numeric.group(2)) if numeric.group(2) else None
    else:
        month_names = {
            "jan": 1, "january": 1, "feb": 2, "february": 2, "mar": 3, "march": 3,
            "apr": 4, "april": 4, "may": 5, "jun": 6, "june": 6, "jul": 7,
            "july": 7, "aug": 8, "august": 8, "sep": 9, "sept": 9, "september": 9,
            "oct": 10, "october": 10, "nov": 11, "november": 11, "dec": 12, "december": 12,
        }
        month_match = re.search(r"(" + "|".join(month_names) + r")", value, re.I)
        if month_match:
            month = month_names[month_match.group(1).casefold()]
            remainder = value.replace(year_match.group(0), " ", 1)
            remainder = re.sub(month_match.group(0), " ", remainder, count=1, flags=re.I)
            day_match = re.search(r"\b(\d{1,2})\b", remainder)
            if day_match and 1 <= int(day_match.group(1)) <= 31:
                day = int(day_match.group(1))

    if month is not None and not 1 <= month <= 12:
        return None
    if day is not None and month is not None:
        import calendar
        if day < 1 or day > calendar.monthrange(year, month)[1]:
            return None
    if year > today.year:
        return "Source reports a publication date after the current year; original metadata retained pending publisher verification."
    if year < today.year:
        return None
    if month is not None and 1 <= month <= 12 and month > today.month:
        return "Source reports a publication date later than today; original metadata retained pending publisher verification."
    if month == today.month and day is not None and day > today.day:
        return "Source reports a publication date later than today; original metadata retained pending publisher verification."
    return None


def _apply_date_validation(article: LiteratureArticle, today: Optional[datetime] = None) -> bool:
    """Set date-validation metadata; return True when a future date needs review."""
    note = _future_date_note(article.publication_date, today=today)
    if note:
        article.metadata["date_validation"] = "FUTURE_DATE_REVIEW"
        article.metadata["date_validation_note"] = note
        return True
    article.metadata.setdefault("date_validation", "NOT_FLAGGED")
    article.metadata.pop("date_validation_note", None)
    return False



def search_pubmed(query: str, max_results: int = 30, date_from: Optional[str] = None,

                  date_to: Optional[str] = None, email: Optional[str] = None) -> list[LiteratureArticle]:

    """Search PubMed and retrieve article metadata plus abstracts via E-utilities."""

    if not query or not query.strip():

        raise ValueError("Search query cannot be empty.")

    if not 1 <= max_results <= MAX_RESULTS_LIMIT:

        raise ValueError(f"max_results must be between 1 and {MAX_RESULTS_LIMIT}.")

    common = {"db": "pubmed", "tool": "Synapse", "email": email or os.getenv("SYNAPSE_NCBI_EMAIL"),
              "api_key": _ncbi_api_key()}

    params = {**common, "term": query.strip(), "retmax": max_results, "retmode": "json", "sort": "relevance"}

    # E-utilities expect YYYY, YYYY/MM or YYYY/MM/DD.
    if date_from:

        params.update({"mindate": date_from.replace("-", "/"), "datetype": "pdat"})

    if date_to:

        params.update({"maxdate": date_to.replace("-", "/"), "datetype": "pdat"})

    search_data = _request_json(_url("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params))

    ids = search_data.get("esearchresult", {}).get("idlist", [])

    if not ids:

        return []

    fetch_params = {**common, "id": ",".join(ids), "rettype": "abstract", "retmode": "xml"}

    xml_text = _request_text(_url("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", fetch_params))

    try:

        root = ET.fromstring(xml_text)

    except ET.ParseError as exc:

        raise LiteratureSourceError("PubMed returned malformed XML.") from exc

    articles: list[LiteratureArticle] = []

    for record in root.findall(".//PubmedArticle"):

        pmid = _clean_text(record.findtext(".//MedlineCitation/PMID"))

        if not pmid:

            continue

        title_node = record.find(".//Article/ArticleTitle")

        title = "".join(title_node.itertext()).strip() if title_node is not None else ""

        authors: list[str] = []

        for author in record.findall(".//Article/AuthorList/Author"):

            collective = _clean_text(author.findtext("CollectiveName"))

            fore = _clean_text(author.findtext("ForeName"))

            last = _clean_text(author.findtext("LastName"))

            name = collective or " ".join(x for x in (fore, last) if x)

            if name:

                authors.append(name)

        abstract_parts: list[str] = []

        for node in record.findall(".//Article/Abstract/AbstractText"):

            text = "".join(node.itertext()).strip()

            label = node.attrib.get("Label") or node.attrib.get("NlmCategory")

            if text:

                abstract_parts.append(f"{label}: {text}" if label and label.upper() != "UNASSIGNED" else text)

        doi = None

        pmcid = None

        for article_id in record.findall(".//PubmedData/ArticleIdList/ArticleId"):

            kind = article_id.attrib.get("IdType", "").lower()

            value = _clean_text(article_id.text)

            if kind == "doi":

                doi = value or None

            elif kind == "pmc":

                pmcid = value or None

        journal = _clean_text(record.findtext(".//Article/Journal/Title")) or None

        pubdate = record.find(".//Article/Journal/JournalIssue/PubDate")

        date_parts: list[str] = []

        if pubdate is not None:

            for tag in ("Year", "Month", "Day", "MedlineDate"):

                value = _clean_text(pubdate.findtext(tag))

                if value:

                    date_parts.append(value)

                    if tag == "MedlineDate" or tag == "Day":

                        break

        publication_date = " ".join(date_parts) or None

        pub_types = [_clean_text(n) for n in record.findall(".//Article/PublicationTypeList/PublicationType") if _clean_text(n)]

        mesh = [_clean_text(n.findtext("DescriptorName")) for n in record.findall(".//MeshHeadingList/MeshHeading") if _clean_text(n.findtext("DescriptorName"))]

        articles.append(LiteratureArticle(

            title=title or f"PubMed article {pmid}", article_id=f"pubmed:{pmid}", doi=doi, pmid=pmid, pmcid=pmcid,

            authors=authors, journal=journal, publication_date=publication_date,

            abstract=_clean_text("\n".join(abstract_parts)) or None, keywords=mesh, publication_types=pub_types,

            source_ids={"PubMed": pmid}, source_urls=[f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"],

            metadata={"source": "PubMed", "abstract_available": bool(abstract_parts)}))

    return articles

def search_europe_pmc(query: str, max_results: int = 30, date_from: Optional[str] = None,

                      date_to: Optional[str] = None, open_access_only: bool = False) -> list[LiteratureArticle]:

    """Search Europe PMC REST API for biomedical records and abstracts."""

    parts = [f"({query.strip()})"]

    if date_from:

        parts.append(f"FIRST_PDATE:[{_date_bound(date_from, upper=False)} TO *]")

    if date_to:

        parts.append(f"FIRST_PDATE:[* TO {_date_bound(date_to, upper=True)}]")

    if open_access_only:

        parts.append("OPEN_ACCESS:Y")

    params = {"query": " AND ".join(parts), "format": "json", "pageSize": max_results,

              "resultType": "core", "sortBy": "RELEVANCE"}

    data = _request_json(_url("https://www.ebi.ac.uk/europepmc/webservices/rest/search", params))

    articles: list[LiteratureArticle] = []

    for item in data.get("resultList", {}).get("result", []):

        title = _clean_text(item.get("title"))

        pmid = _clean_text(item.get("pmid")) or None

        pmcid = _clean_text(item.get("pmcid")) or None

        doi = _clean_text(item.get("doi")) or None

        author_string = _clean_text(item.get("authorString"))

        authors = [a.strip() for a in author_string.split(",") if a.strip()] if author_string else []

        identifier = pmid or pmcid or doi or _clean_text(item.get("id"))

        if not title or not identifier:

            continue

        url = (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else

               f"https://europepmc.org/article/MED/{pmid}" if pmid else

               f"https://europepmc.org/article/PMC/{pmcid}" if pmcid else

               f"https://doi.org/{doi}" if doi else "https://europepmc.org/")

        articles.append(LiteratureArticle(

            title=title, article_id=f"europepmc:{identifier}", doi=doi, pmid=pmid, pmcid=pmcid,

            authors=authors, journal=_clean_text(item.get("journalTitle")) or None,

            publication_date=_clean_text(item.get("firstPublicationDate") or item.get("pubYear")) or None,

            abstract=_clean_text(item.get("abstractText")) or None,

            keywords=[_clean_text(x) for x in item.get("keywordList", {}).get("keyword", []) if _clean_text(x)],

            source_ids={"Europe PMC": identifier}, source_urls=[url],

            is_open_access=(str(item.get("isOpenAccess", "")).upper() == "Y"),

            citation_count=int(item["citedByCount"]) if str(item.get("citedByCount", "")).isdigit() else None,

            metadata={"source": "Europe PMC", "abstract_available": bool(item.get("abstractText"))}))

    return articles

def search_openalex(query: str, max_results: int = 30, date_from: Optional[str] = None,

                    date_to: Optional[str] = None, email: Optional[str] = None) -> list[LiteratureArticle]:

    """Search OpenAlex works. Abstracts are reconstructed from its inverted index when available."""

    filters: list[str] = []

    # OpenAlex needs full YYYY-MM-DD dates. Previously "YYYY-MM" became the
    # invalid "YYYY-MM-01-01" (audit finding); _date_bound expands correctly.
    if date_from:

        filters.append(f"from_publication_date:{_date_bound(date_from, upper=False)}")

    if date_to:

        filters.append(f"to_publication_date:{_date_bound(date_to, upper=True)}")

    params: dict[str, Any] = {"search": query.strip(), "per-page": max_results, "sort": "relevance_score:desc"}

    if filters:

        params["filter"] = ",".join(filters)

    if email or os.getenv("SYNAPSE_OPENALEX_EMAIL"):

        params["mailto"] = email or os.getenv("SYNAPSE_OPENALEX_EMAIL")

    data = _request_json(_url("https://api.openalex.org/works", params))

    articles: list[LiteratureArticle] = []

    for item in data.get("results", []):

        title = _clean_text(item.get("title"))

        work_id = _clean_text(item.get("id"))

        if not title or not work_id:

            continue

        doi = _clean_text(item.get("doi")) or None

        doi = re.sub(r"^https?://doi\.org/", "", doi, flags=re.I) if doi else None
        authors = [_clean_text(a.get("author", {}).get("display_name")) for a in item.get("authorships", [])]

        authors = [a for a in authors if a]

        inverted = item.get("abstract_inverted_index") or {}

        words: list[tuple[int, str]] = []

        for word, positions in inverted.items():

            if isinstance(positions, list):

                words.extend((int(pos), word) for pos in positions if isinstance(pos, int))

        abstract = " ".join(word for _, word in sorted(words)) or None

        source_ids = {"OpenAlex": work_id}

        ids = item.get("ids") or {}

        pmid_url = ids.get("pmid")

        pmid_match = re.search(r"(\d+)(?:/?$)", pmid_url or "")

        pmid = pmid_match.group(1) if pmid_match else None

        if pmid:

            source_ids["PubMed"] = pmid

        venue = item.get("primary_location") or {}

        source = venue.get("source") or {}

        journal = _clean_text(source.get("display_name")) or None

        date = _clean_text(item.get("publication_date")) or None

        cited = item.get("cited_by_count")

        oa = item.get("open_access") or {}

        url = _clean_text(item.get("doi")) or work_id

        articles.append(LiteratureArticle(

            title=title, article_id=f"openalex:{work_id.rsplit('/', 1)[-1]}", doi=doi, pmid=pmid,

            authors=authors, journal=journal, publication_date=date, abstract=abstract,

            source_ids=source_ids, source_urls=[url], is_open_access=oa.get("is_oa") if isinstance(oa.get("is_oa"), bool) else None,

            citation_count=int(cited) if isinstance(cited, int) else None,

            metadata={"source": "OpenAlex", "openalex_id": work_id, "abstract_available": bool(abstract)}))

    return articles

def search_crossref(query: str, max_results: int = 30, date_from: Optional[str] = None,

                    date_to: Optional[str] = None, email: Optional[str] = None) -> list[LiteratureArticle]:

    """Search Crossref works for DOI and bibliographic metadata."""

    filters: list[str] = []

    if date_from:

        filters.append("from-pub-date:" + date_from)

    if date_to:

        filters.append("until-pub-date:" + date_to)

    params: dict[str, Any] = {"query": query.strip(), "rows": max_results, "select": "DOI,title,author,container-title,published-print,published-online,published,abstract,type,URL"}

    if filters:

        params["filter"] = ",".join(filters)

    contact = email or os.getenv("SYNAPSE_CROSSREF_EMAIL")

    if contact:

        params["mailto"] = contact

    data = _request_json(_url("https://api.crossref.org/works", params))

    articles: list[LiteratureArticle] = []

    for item in data.get("message", {}).get("items", []):

        title = _clean_text((item.get("title") or [""])[0])

        doi = _clean_text(item.get("DOI")) or None

        if not title:

            continue

        authors = []

        for author in item.get("author", []):

            name = " ".join(x for x in (_clean_text(author.get("given")), _clean_text(author.get("family"))) if x)

            if name:

                authors.append(name)

        date_parts = item.get("published-print") or item.get("published-online") or item.get("published") or {}

        date_values = date_parts.get("date-parts", [[]])

        date = "-".join(str(x).zfill(2) if i else str(x) for i, x in enumerate(date_values[0])) if date_values and date_values[0] else None

        abstract = _clean_text(item.get("abstract")) or None

        articles.append(LiteratureArticle(

            title=title, article_id=f"crossref:{doi or title}", doi=doi, authors=authors,

            journal=_clean_text((item.get("container-title") or [""])[0]) or None,

            publication_date=date, abstract=abstract, source_ids={"Crossref": doi or title},

            source_urls=[_clean_text(item.get("URL")) or (f"https://doi.org/{doi}" if doi else "https://search.crossref.org/")],

            publication_types=[_clean_text(item.get("type"))] if item.get("type") else [],

            metadata={"source": "Crossref", "abstract_available": bool(abstract)}))

    return articles

def _is_placeholder_journal(value: Optional[str]) -> bool:
    """Return True when a source/platform label was mistakenly stored as a journal."""
    if not value:
        return False
    normalized = re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()
    return normalized in {"pubmed", "europe pmc", "openalex", "crossref", "ncbi", "pubmed central"}

def deduplicate_articles(articles: list[LiteratureArticle]) -> list[LiteratureArticle]:
    """Merge records conservatively: stable identifiers outrank title-only matches.

    Different DOI/PMID/PMCID records are never merged solely because their titles
    happen to match. Title fallback is used only when neither record has a stable ID.
    """
    merged: list[LiteratureArticle] = []
    strong_key_to_index: dict[str, int] = {}
    title_only_to_index: dict[str, int] = {}

    def strong_keys(article: LiteratureArticle) -> list[str]:
        return [key for key in article.identity_keys() if key.startswith(("doi:", "pmid:", "pmcid:"))]

    def title_key(article: LiteratureArticle) -> Optional[str]:
        normalized = re.sub(r"[^a-z0-9]+", " ", (article.title or "").lower()).strip()
        return "title:" + normalized if normalized else None

    for article in articles:
        ids = strong_keys(article)
        tkey = title_key(article)
        matched_index = next((strong_key_to_index[key] for key in ids if key in strong_key_to_index), None)
        if matched_index is None and not ids and tkey and tkey in title_only_to_index:
            matched_index = title_only_to_index[tkey]

        if matched_index is None:
            merged.append(article)
            matched_index = len(merged) - 1
        else:
            target = merged[matched_index]
            # Fill missing bibliographic fields, never overwrite existing identifiers.
            for attr in ("doi", "pmid", "pmcid", "journal", "publication_date", "citation_count", "is_open_access"):
                incoming = getattr(article, attr)
                current = getattr(target, attr)
                if attr == "journal":
                    if _is_placeholder_journal(current):
                        setattr(target, attr, None)
                        current = None
                    if _is_placeholder_journal(incoming):
                        incoming = None
                if current in (None, "", 0) and incoming not in (None, "", 0):
                    setattr(target, attr, incoming)
            if len(article.abstract or "") > len(target.abstract or ""):
                target.abstract = article.abstract
            if article.title and (not target.title or target.title.lower().startswith("untitled")):
                target.title = article.title
            target.authors = _deduplicate_author_names(target.authors + article.authors)
            target.keywords = _normalize_keywords(target.keywords + article.keywords)
            target.publication_types = list(dict.fromkeys(target.publication_types + article.publication_types))
            if len(article.title or "") > len(target.title or ""):
                target.title = article.title
            if article.journal and not _is_placeholder_journal(article.journal):
                if not target.journal or _is_placeholder_journal(target.journal) or len(article.journal) > len(target.journal):
                    target.journal = article.journal
            target.source_ids.update(article.source_ids)
            target.source_urls = list(dict.fromkeys(target.source_urls + article.source_urls))
            target.metadata.update({k: v for k, v in article.metadata.items() if k not in target.metadata or not target.metadata[k]})

        target = merged[matched_index]
        target.abstract = _clean_text(target.abstract) or None if target.abstract else None
        target.title = _clean_text(target.title)
        target.journal = _clean_text(target.journal) or None if target.journal else None
        if _is_placeholder_journal(target.journal):
            target.journal = None
        target.authors = _deduplicate_author_names(target.authors)
        target.keywords = _normalize_keywords(target.keywords)
        for key in strong_keys(target):
            strong_key_to_index[key] = matched_index
        # Only index title-only records without stable identifiers. This avoids
        # collapsing distinct publications that share a generic title.
        if not strong_keys(target) and tkey:
            title_only_to_index[tkey] = matched_index

    return merged

def _tokenize(text: str) -> list[str]:
    """Tokenize biomedical text and ignore Boolean query operators."""
    stop = {
        "the", "and", "for", "with", "from", "that", "this", "into", "using",
        "study", "based", "are", "was", "were", "its", "their", "about", "via",
        "not", "but", "can", "may", "new", "of", "in", "on", "to", "a", "an",
        "by", "as", "at", "or", "we", "our", "these", "those", "which", "how",
        "and/or", "not", "near", "adj", "all", "any",
    }
    tokens = re.findall(r"[a-zA-Z0-9]+(?:[.+-][a-zA-Z0-9]+)*", text.lower())
    return [token for token in tokens if len(token) > 1 and token not in stop]


def _coverage(query_tokens: list[str], text: str) -> float:
    if not query_tokens:
        return 0.0
    text_tokens = set(_tokenize(text))
    unique_query_tokens = set(query_tokens)
    return len(unique_query_tokens & text_tokens) / len(unique_query_tokens)


def _query_phrase(query: str) -> str:
    """Return a conservative phrase candidate without Boolean operators."""
    value = re.sub(r"[()\[\]{}]", " ", query.lower())
    value = re.sub(r"\b(and|or|not)\b", " ", value)
    value = value.replace("\"", "").replace("'", "")
    return re.sub(r"\s+", " ", value).strip()


def rank_articles(
    articles: list[LiteratureArticle],
    query: str,
    biological_context: Optional[dict[str, Any]] = None,
) -> list[LiteratureArticle]:
    """Rank query relevance with transparent field weights and coverage diagnostics.

    This is an information-retrieval score, not a measure of evidence quality,
    methodological rigor, or truth. Citation counts are deliberately not score inputs.
    """
    query_clean = re.sub(r"\s+", " ", query).strip().lower()
    query_tokens = _tokenize(query_clean)
    query_token_set = set(query_tokens)
    context_values: list[str] = []
    for value in (biological_context or {}).values():
        if isinstance(value, (list, tuple, set)):
            context_values.extend(str(item) for item in value if item)
        elif value:
            context_values.append(str(value))
    context_tokens = set(_tokenize(" ".join(context_values))) - query_token_set

    if not query_token_set:
        for article in articles:
            article.relevance_score = 0.0
            article.relevance_components = {
                "title_coverage": 0.0, "abstract_coverage": 0.0,
                "keyword_coverage": 0.0, "query_term_coverage": 0.0,
                "exact_phrase_bonus": 0.0,
            }
        return articles

    phrase = _query_phrase(query)
    for article in articles:
        title = article.title or ""
        abstract = article.abstract or ""
        keywords = " ".join(article.keywords)
        title_cov = _coverage(query_tokens, title)
        abstract_cov = _coverage(query_tokens, abstract)
        keyword_cov = _coverage(query_tokens, keywords)
        all_text = f"{title} {abstract} {keywords}"
        all_tokens = set(_tokenize(all_text))
        term_coverage = len(query_token_set & all_tokens) / len(query_token_set)
        title_normalized = re.sub(r"[^a-z0-9]+", " ", title.lower()).strip()
        phrase_normalized = re.sub(r"[^a-z0-9]+", " ", phrase).strip()
        exact_phrase = bool(phrase_normalized and phrase_normalized in title_normalized)
        context_text = all_text
        context_cov = _coverage(sorted(context_tokens), context_text) if context_tokens else 0.0
        # Field weights sum to 100. Context is a small, separately disclosed bonus.
        base_score = 50.0 * title_cov + 35.0 * abstract_cov + 15.0 * keyword_cov
        phrase_bonus = 5.0 if exact_phrase else 0.0
        context_bonus = 5.0 * context_cov
        score = min(100.0, base_score + phrase_bonus + context_bonus)
        matched_terms = sorted(query_token_set & all_tokens)
        missing_terms = sorted(query_token_set - all_tokens)
        article.relevance_score = round(score, 2)
        article.relevance_components = {
            "title_coverage": round(title_cov, 4),
            "abstract_coverage": round(abstract_cov, 4),
            "keyword_coverage": round(keyword_cov, 4),
            "query_term_coverage": round(term_coverage, 4),
            "exact_phrase_bonus": phrase_bonus,
            "biological_context_coverage": round(context_cov, 4),
            "matched_query_terms": matched_terms,
            "missing_query_terms": missing_terms,
            "matched_context_terms": sorted(context_tokens & all_tokens),
            "ranking_method": "title 50%, abstract 35%, keywords 15%, phrase bonus up to 5, context bonus up to 5; capped at 100",
        }
        # Platform labels are not journal titles and must not inflate metadata completeness.
        if _is_placeholder_journal(article.journal):
            article.journal = None
        # Data completeness is shown separately and never masquerades as evidence quality.
        present = sum(bool(value) for value in (article.title, article.authors, article.journal,
                       article.publication_date, article.abstract, article.doi or article.pmid,
                       article.keywords))
        article.metadata["metadata_completeness"] = round(present / 7.0, 3)
        missing = []
        for field_name, value in (("authors", article.authors), ("journal", article.journal),
                                  ("publication_date", article.publication_date),
                                  ("abstract", article.abstract), ("DOI/PMID", article.doi or article.pmid)):
            if not value:
                missing.append(field_name)
        article.metadata["missing_core_metadata"] = missing

    return sorted(articles, key=lambda article: (
        article.relevance_score if article.relevance_score is not None else -1,
        article.title.lower(),
    ), reverse=True)


def _journal_from_pubmed_xml(xml_text: str) -> Optional[str]:
    """Extract the canonical journal title from a PubMed efetch XML response."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return None
    for record in root.findall(".//PubmedArticle"):
        title = _clean_text(record.findtext(".//Article/Journal/Title"))
        if title and not _is_placeholder_journal(title):
            return title
    return None


def _publication_types_from_pubmed_xml(xml_text: str) -> list[str]:
    """Extract normalized NLM publication types from a PubMed XML response."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []
    values: list[str] = []
    for node in root.findall(".//Article/PublicationTypeList/PublicationType"):
        value = _clean_text("".join(node.itertext()))
        if value and value.casefold() not in {item.casefold() for item in values}:
            values.append(value)
    return values


def _enrich_missing_metadata(articles: list[LiteratureArticle], limit: int = 12) -> tuple[int, list[str]]:
    """Fill missing bibliographic fields from Europe PMC and Crossref.

    Only records with a DOI or PMID and at least one missing field are queried.
    Enrichment is deliberately bounded and failures never discard search results.
    """
    enriched = 0
    warnings: list[str] = []
    candidates = [
        article for article in articles
        if (article.pmid or article.normalized_doi())
        and (
            not article.abstract
            or not article.journal
            or not article.keywords
            or not article.authors
            or not article.publication_date
            or not article.publication_types
        )
    ][:max(0, limit)]

    for article in candidates:
        changed = False
        enrichment_sources: list[str] = []
        query_terms: list[str] = []
        if article.pmid:
            query_terms.append(f"EXT_ID:{article.pmid}")
        elif article.normalized_doi():
            query_terms.append(f'DOI:"{article.normalized_doi()}"')

        if query_terms:
            try:
                params = {
                    "query": " OR ".join(query_terms),
                    "format": "json",
                    "pageSize": 1,
                    "resultType": "core",
                }
                payload = _request_json(
                    _url("https://www.ebi.ac.uk/europepmc/webservices/rest/search", params),
                    timeout=min(DEFAULT_TIMEOUT, 12),
                )
                results = payload.get("resultList", {}).get("result", [])
                if results:
                    record = results[0]
                    enrichment_sources.append("Europe PMC")
                    fields = {
                        "journal": _clean_text(record.get("journalTitle")) or None,
                        "publication_date": _clean_text(
                            record.get("firstPublicationDate") or record.get("pubYear")
                        ) or None,
                        "abstract": _clean_text(record.get("abstractText")) or None,
                    }
                    for field_name, value in fields.items():
                        if value and not getattr(article, field_name):
                            setattr(article, field_name, value)
                            changed = True

                    keyword_values = record.get("keywordList", {}).get("keyword", [])
                    if not article.keywords and isinstance(keyword_values, list):
                        article.keywords = list(dict.fromkeys(
                            _clean_text(item) for item in keyword_values if _clean_text(item)
                        ))
                        changed = changed or bool(article.keywords)

                    author_string = _clean_text(record.get("authorString"))
                    if not article.authors and author_string:
                        article.authors = [name.strip() for name in author_string.split(",") if name.strip()]
                        changed = changed or bool(article.authors)

                    if record.get("pmcid") and not article.pmcid:
                        article.pmcid = _clean_text(record["pmcid"]) or None
                        changed = changed or bool(article.pmcid)
                    if str(record.get("isOpenAccess", "")).upper() == "Y" and article.is_open_access is None:
                        article.is_open_access = True
                        changed = True
            except (LiteratureSourceError, ValueError, TypeError, KeyError) as exc:
                if _http_error_is_not_found(exc):
                    warnings.append(f"Europe PMC has no matching metadata record for '{article.title[:70]}'.")
                else:
                    warnings.append(f"Metadata enrichment via Europe PMC failed for '{article.title[:70]}': {exc}")

        # PubMed efetch is a reliable final fallback for canonical journal titles
        # when a merged record had a repository name mistakenly stored as journal.
        if article.pmid and (not article.journal or not article.publication_types):
            try:
                params = {
                    "db": "pubmed",
                    "id": article.pmid,
                    "rettype": "abstract",
                    "retmode": "xml",
                    "tool": "Synapse",
                    "email": os.getenv("SYNAPSE_NCBI_EMAIL"),
                }
                xml_text = _request_text(
                    _url("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params),
                    timeout=min(DEFAULT_TIMEOUT, 12),
                )
                journal_title = _journal_from_pubmed_xml(xml_text)
                if not article.journal and journal_title:
                    article.journal = journal_title
                    enrichment_sources.append("PubMed")
                    changed = True
                pub_types = _publication_types_from_pubmed_xml(xml_text)
                if not article.publication_types and pub_types:
                    article.publication_types = pub_types
                    enrichment_sources.append("PubMed")
                    changed = True
            except (LiteratureSourceError, ValueError, TypeError, KeyError) as exc:
                warnings.append(
                    f"PubMed metadata lookup failed for PMID {article.pmid}: {exc}"
                )

        # Crossref is a secondary fallback for DOI-bearing records, particularly
        # for journal and publication metadata. Its abstracts may contain JATS XML.
        if article.normalized_doi() and (
            not article.abstract or not article.journal or not article.publication_date
        ):
            try:
                doi = urllib.parse.quote(article.normalized_doi(), safe="/")
                payload = _request_json(
                    f"https://api.crossref.org/works/{doi}",
                    timeout=min(DEFAULT_TIMEOUT, 12),
                    headers={"Accept": "application/json"},
                )
                record = payload.get("message", {})
                enrichment_sources.append("Crossref")
                if not article.journal:
                    container = record.get("container-title") or []
                    if container:
                        article.journal = _clean_text(container[0]) or None
                        changed = changed or bool(article.journal)
                if not article.publication_date:
                    date_parts = (
                        record.get("published-print", {}).get("date-parts")
                        or record.get("published-online", {}).get("date-parts")
                        or record.get("created", {}).get("date-parts")
                        or []
                    )
                    if date_parts and date_parts[0]:
                        article.publication_date = "-".join(str(x) for x in date_parts[0][:3])
                        changed = True
                if not article.abstract and record.get("abstract"):
                    article.abstract = _clean_text(record["abstract"]) or None
                    changed = changed or bool(article.abstract)
                if not article.authors and record.get("author"):
                    names = []
                    for author in record["author"]:
                        name = " ".join(
                            part for part in (
                                _clean_text(author.get("given")),
                                _clean_text(author.get("family")),
                            ) if part
                        )
                        if name:
                            names.append(name)
                    article.authors = names
                    changed = changed or bool(names)
            except (LiteratureSourceError, ValueError, TypeError, KeyError) as exc:
                if _http_error_is_not_found(exc):
                    warnings.append(f"Crossref has no matching DOI record for '{article.title[:70]}'; existing metadata was retained.")
                else:
                    warnings.append(f"Metadata enrichment via Crossref failed for '{article.title[:70]}': {exc}")

        # Normalize repeated labels after cross-source metadata merging.
        if article.authors:
            normalized_authors = _deduplicate_author_names(article.authors)
            if normalized_authors != article.authors:
                article.authors = normalized_authors
                changed = True
        if article.keywords:
            normalized_keywords = _normalize_keywords(article.keywords)
            if normalized_keywords != article.keywords:
                article.keywords = normalized_keywords
                changed = True

        if changed:
            article.metadata["metadata_enriched"] = True
            article.metadata["metadata_enrichment_sources"] = list(dict.fromkeys(enrichment_sources))
            enriched += 1

    return enriched, warnings


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

def article_to_dict(article: LiteratureArticle) -> dict[str, Any]:

    return asdict(article)

def export_results(run: LiteratureSearchRun, path: str | Path) -> Path:

    """Export a complete search run as JSON, CSV, TXT, RIS, or BibTeX."""

    destination = Path(path).expanduser()

    suffix = destination.suffix.lower()

    if suffix not in {".json", ".csv", ".txt", ".ris", ".bib", ".bibtex"}:

        raise ValueError("Supported export formats: .json, .csv, .txt, .ris, .bib")

    destination.parent.mkdir(parents=True, exist_ok=True)

    if suffix == ".json":

        payload = {"module": MODULE_NAME, "module_version": MODULE_VERSION, "run": {

            "run_id": run.run_id, "started_at_utc": run.started_at_utc, "completed_at_utc": run.completed_at_utc,

            "query": asdict(run.query), "source_status": run.source_status, "warnings": run.warnings,

            "errors": run.errors, "total_results_before_deduplication": run.total_results_before_deduplication,

            "total_results_after_deduplication": run.total_results_after_deduplication,

            "articles": [article_to_dict(a) for a in run.articles],

            "ranking_method": "Transparent weighted text overlap: title 50%, abstract 35%, keywords 15%, phrase bonus up to 5 points, optional biological-context bonus up to 5 points; capped at 100. Citation counts are not score inputs.",

            "scientific_note": "Relevance scores estimate query matching, not evidence quality or truth."}}

        destination.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    elif suffix == ".csv":

        fields = ["title", "doi", "pmid", "pmcid", "authors", "journal", "publication_date", "abstract", "keywords", "publication_types", "source_ids", "source_urls", "is_open_access", "citation_count", "relevance_score", "relevance_components", "metadata"]

        with destination.open("w", newline="", encoding="utf-8-sig") as handle:

            writer = csv.DictWriter(handle, fieldnames=fields)

            writer.writeheader()

            for article in run.articles:

                row = article_to_dict(article)

                for key in ("authors", "keywords", "publication_types", "source_ids", "source_urls", "relevance_components", "metadata"):

                    row[key] = json.dumps(row[key], ensure_ascii=False)

                writer.writerow({k: row.get(k) for k in fields})

    elif suffix == ".txt":

        lines = [MODULE_NAME, f"Run ID: {run.run_id}", f"Query: {run.query.query}", f"Sources: {', '.join(run.query.sources)}", f"Unique articles: {len(run.articles)}", "", "SCIENTIFIC NOTE: Relevance scores are not evidence-quality scores.", ""]

        for i, a in enumerate(run.articles, 1):

            lines += [f"{i}. {a.title}", f"   Authors: {', '.join(a.authors) or 'Not available'}", f"   Journal/date: {a.journal or 'Unknown'} / {a.publication_date or 'Unknown'}" + (" [REVIEW: FUTURE DATE]" if a.metadata.get("date_validation") == "FUTURE_DATE_REVIEW" else ""), f"   DOI: {a.doi or 'Not available'} | PMID: {a.pmid or 'Not available'}", f"   Relevance: {a.relevance_score if a.relevance_score is not None else 'N/A'}", f"   URL: {a.source_urls[0] if a.source_urls else 'Not available'}", f"   Abstract: {a.abstract or 'Not available'}", ""]

        destination.write_text("\n".join(lines), encoding="utf-8")

    elif suffix == ".ris":

        lines = []

        for a in run.articles:

            lines.append("TY  - JOUR")

            lines.append("TI  - " + a.title)

            for author in a.authors:

                lines.append("AU  - " + author)

            if a.journal: lines.append("JO  - " + a.journal)

            if a.publication_date: lines.append("PY  - " + a.publication_date[:4])

            if a.doi: lines.append("DO  - " + a.doi)

            if a.pmid: lines.append("AN  - PMID:" + a.pmid)

            if a.abstract: lines.append("AB  - " + a.abstract.replace("\n", " "))
            if a.metadata.get("date_validation") == "FUTURE_DATE_REVIEW": lines.append("N1  - REVIEW: FUTURE PUBLICATION DATE. " + a.metadata.get("date_validation_note", "Verify source date."))

            for url in a.source_urls: lines.append("UR  - " + url)

            lines.append("ER  - ")

            lines.append("")

        destination.write_text("\n".join(lines), encoding="utf-8")

    else:

        entries = []

        for index, a in enumerate(run.articles, 1):

            key = re.sub(r"[^A-Za-z0-9]", "", (a.authors[0].split()[-1] if a.authors else "article") + (a.publication_date or "")[:4] + str(index)) or f"article{index}"

            def bib_escape(value: str) -> str:
                escaped: list[str] = []
                for character in value:
                    if character == "\\":
                        escaped.append(r"\textbackslash{}")
                    elif character == "{":
                        escaped.append(r"\{")
                    elif character == "}":
                        escaped.append(r"\}")
                    else:
                        escaped.append(character)
                return "".join(escaped)

            fields = [f"  title = {{{bib_escape(a.title)}}}"]

            if a.authors: fields.append("  author = {" + " and ".join(bib_escape(x) for x in a.authors) + "}")

            if a.journal: fields.append(f"  journal = {{{bib_escape(a.journal)}}}")

            if a.publication_date: fields.append(f"  year = {{{a.publication_date[:4]}}}")

            if a.doi: fields.append(f"  doi = {{{bib_escape(a.doi)}}}")

            if a.abstract: fields.append(f"  abstract = {{{bib_escape(a.abstract)}}}")
            if a.metadata.get("date_validation") == "FUTURE_DATE_REVIEW": fields.append("  note = {" + bib_escape("REVIEW: FUTURE PUBLICATION DATE. " + a.metadata.get("date_validation_note", "Verify source date.")) + "}")

            entries.append("@article{" + key + ",\n" + ",\n".join(fields) + "\n}")

        destination.write_text("\n\n".join(entries), encoding="utf-8")

    return destination

def print_results(run: LiteratureSearchRun, limit: int = 10) -> None:

    style.header('Literature search results')

    style.field('Query', f"{run.query.query}")

    style.field('Sources', f"{', '.join(run.query.sources)}")

    print(f"Retrieved before deduplication: {run.total_results_before_deduplication}")

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

        print("WARNINGS")

        for warning in run.warnings: print(f"- {warning}")

    if run.errors:

        print("SOURCE ERRORS")

        for error in run.errors: print(f"- {error}")

    print("\nInterpretation: relevance ranking reflects text/query overlap, not scientific quality or truth.")

def print_article_detail(run: LiteratureSearchRun, article_index: int) -> None:

    """Display complete normalized metadata for one result (1-based index)."""

    if article_index < 1 or article_index > len(run.articles):

        print("Article number is outside the result range.")

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

    print(f"DOI: {article.doi or 'Not available'}")

    print(f"PMID: {article.pmid or 'Not available'} | PMCID: {article.pmcid or 'Not available'}")

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

    print("\nABSTRACT\n--------")

    print(article.abstract or "Abstract not available from the retrieved sources.")

    print("\nKEYWORDS / SUBJECT TERMS")

    print(", ".join(_normalize_keywords(article.keywords)) or "Not available")

    print("\nSCIENTIFIC LIMITATION: metadata and relevance are not an appraisal of study quality or validity.")

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

    check("relevance_ranking", ranked[0].relevance_score >= ranked[1].relevance_score)

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

    check("biological_context_match_recorded", "zebrafish" in context_ranked[0].relevance_components["matched_context_terms"])

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

    style.field('Self-test', f"{sum(ok for _, ok in tests)}/{len(tests)} passed")

    return all(ok for _, ok in tests)

def build_sequence_literature_context(header: Optional[str] = None, orf: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    """Derive an optional search suggestion from the active sequence.

    Only the free-text sequence header (and the ORF label, if any) is used;
    no sequence data leaves the program. The suggestion is shown to the user
    and is never submitted without their confirmation.
    """
    context: dict[str, Any] = {}
    words = [w for w in re.split(r"[\s|,;_]+", (header or "").strip()) if w]
    words = [w for w in words if not re.fullmatch(r"(?i)(topology|len|length)=.*|\d+", w)]
    if words:
        context["sequence_header"] = header.strip()
        context["suggested_query"] = " ".join(words[:6])
    if orf:
        context["orf"] = f"rank {orf.get('rank', '?')} frame {orf.get('frame', '?')}"
    return context


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

            print("Optional environment variables: SYNAPSE_NCBI_EMAIL, SYNAPSE_OPENALEX_EMAIL, SYNAPSE_CROSSREF_EMAIL.")

            continue

        if choice != "1":

            style.invalid_selection()

            continue

        print("Tip: use specific biological terms (gene/protein, organism, disease, pathway) instead of broad single words.")

        print("Example: BRCA1 DNA repair breast cancer Homo sapiens")

        if suggested:
            print(f"Suggestion from the active sequence header: {suggested}")
            search_text = input("Search query [Enter = use suggestion]: ").strip() or suggested
        else:
            search_text = input("Search query (Boolean terms supported by each source): ").strip()

        if not search_text:

            print("Query cannot be empty.")

            continue

        max_raw = style.ask('Maximum results per source', default='30')

        try:

            max_results = int(max_raw) if max_raw else 30

            source_raw = input("Sources [1=PubMed, 2=Europe PMC, 3=OpenAlex, 4=Crossref, Enter=1,2,3]: ").strip()

            source_map = {"1": "pubmed", "2": "europe_pmc", "3": "openalex", "4": "crossref"}

            if source_raw:

                selected = [source_map[x.strip()] for x in source_raw.split(",") if x.strip() in source_map]

            else:

                selected = ["pubmed", "europe_pmc", "openalex"]

            if not selected:

                print("Select at least one valid source.")

                continue

            date_from = style.ask('From date YYYY or YYYY-MM-DD', default='none') or None

            date_to = style.ask('To date YYYY or YYYY-MM-DD', default='none') or None

            oa = input("Open-access only? (y/N): ").strip().lower() == "y"

            sort = input("Sort [relevance/date/citations] (relevance): ").strip().lower() or "relevance"

            query = LiteratureQuery(search_text, selected, max_results, date_from, date_to, oa, sort,
                                    biological_context={k: v for k, v in context.items() if k == "sequence_header"})

            query.validate()

        except (ValueError, KeyError) as exc:

            style.field('Invalid configuration', f"{exc}")

            continue

        print("\nSearching selected literature sources. A source may fail independently...\n")

        result = run_literature_search(query)

        if not result.run:

            print(result.message or "Search did not produce a run.")

            continue

        run = result.run

        if sort == "date":

            run.articles.sort(key=lambda a: a.publication_date or "", reverse=True)

        elif sort == "citations":

            run.articles.sort(key=lambda a: a.citation_count or 0, reverse=True)

        print_results(run)

        while run.articles:

            detail_choice = input("View full article detail? Enter result number, or press Enter to continue: ").strip()

            if not detail_choice:

                break

            try:

                print_article_detail(run, int(detail_choice))

            except ValueError:

                print("Enter a valid result number.")

        export = input("Export results? Enter path with .json/.csv/.txt/.ris/.bib, or press Enter to skip: ").strip()

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

            style.field('Search failed', f"{exc}")

            return 1

    _interactive_menu()

    return 0

if __name__ == "__main__":

    raise SystemExit(main())
