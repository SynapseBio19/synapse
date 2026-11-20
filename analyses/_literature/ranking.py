"""De-duplication of records and the transparent relevance ranking."""

from __future__ import annotations

import re
from typing import Any, Optional

from analyses._literature.models import LiteratureArticle
from analyses._literature.normalize import (
    _clean_text,
    _deduplicate_author_names,
    _normalize_keywords,
)


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
