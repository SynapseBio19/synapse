"""Text, person-name, keyword and date normalization and checks."""

from __future__ import annotations

import html
import re
import unicodedata
from datetime import datetime, timezone
from typing import Any, Optional

from analyses._literature.models import LiteratureArticle, LiteratureQuery


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
    keys: list[list[tuple[str, tuple[str, ...], str]]] = []

    def key(family: str, given: str) -> tuple[str, tuple[str, ...], str]:
        family_key = _normalize_person_name(family)
        given_tokens = tuple(_normalize_person_name(given).split())
        return family_key, given_tokens, " ".join(given_tokens)

    def readings(name: str) -> list[tuple[str, tuple[str, ...], str]]:
        """The ways to read a name as (family, given tokens, given text).

        "Smith, Ann" and "Ann Smith" have one reading. A name ending in one
        to three capital initials ("Smith A", "Smith AB", "Consortium T.",
        the form of Europe PMC's authorString) is also read as family name
        followed by initials; the first reading stays, because a family
        name written in capitals ("Wei LI") looks the same.
        """
        cleaned = re.sub(r"\s+", " ", str(name or "")).strip(" ,;")
        if not cleaned:
            return [("", (), "")]
        if "," in cleaned:
            family, given = [part.strip() for part in cleaned.split(",", 1)]
            return [key(family, given)]
        tokens = cleaned.split()
        found = [key(tokens[-1], " ".join(tokens[:-1]))]
        if len(tokens) >= 2 and re.fullmatch(r"(?:[A-Z]\.?){1,3}", tokens[-1]):
            found.append(key(" ".join(tokens[:-1]), " ".join(tokens[-1].replace(".", ""))))
        return found

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
        name_readings = readings(name)
        if not name_readings[0][0]:
            continue
        duplicate_index = None
        for index, existing_readings in enumerate(keys):
            for family_key, given_tokens, normalized_given in name_readings:
                for existing_family, existing_given_tokens, existing_normalized in existing_readings:
                    if family_key != existing_family:
                        continue
                    if (normalized_given == existing_normalized
                            or compatible_initials(given_tokens, existing_given_tokens)):
                        duplicate_index = index
                        break
                if duplicate_index is not None:
                    break
            if duplicate_index is not None:
                break

        if duplicate_index is None:
            output.append(name)
            keys.append(name_readings)
        elif len(name) > len(output[duplicate_index]):
            output[duplicate_index] = name
            keys[duplicate_index] = name_readings

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
