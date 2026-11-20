"""Filling missing metadata from Europe PMC, PubMed and Crossref."""

from __future__ import annotations

import os
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Optional

from analyses._literature.models import DEFAULT_TIMEOUT, LiteratureArticle, LiteratureSourceError
from analyses._literature.normalize import (
    _clean_text,
    _deduplicate_author_names,
    _normalize_keywords,
)
from analyses._literature.http import _http_error_is_not_found, _request_json, _request_text, _url
from analyses._literature.ranking import _is_placeholder_journal


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
        normalized_doi = article.normalized_doi()
        if normalized_doi and (
            not article.abstract or not article.journal or not article.publication_date
        ):
            try:
                doi = urllib.parse.quote(normalized_doi, safe="/")
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
