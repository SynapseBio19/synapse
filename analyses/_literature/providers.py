"""The PubMed, Europe PMC, OpenAlex and Crossref searches."""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from typing import Any, Optional

from analyses._literature.models import LiteratureArticle, LiteratureSourceError, MAX_RESULTS_LIMIT
from analyses._literature.normalize import _clean_text, _date_bound
from analyses._literature.http import _ncbi_api_key, _request_json, _request_text, _url


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
