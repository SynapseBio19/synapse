"""JSON, CSV, text, RIS and BibTeX exports."""

from __future__ import annotations

import csv
import json
import re
from core.textformat import field_line
from dataclasses import asdict
from pathlib import Path
from typing import Any

from analyses._literature.models import (
    LiteratureArticle,
    LiteratureSearchRun,
    MODULE_NAME,
    MODULE_VERSION,
)


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

        lines = [MODULE_NAME, field_line('Run ID', f"{run.run_id}"), field_line('Query', f"{run.query.query}"), field_line('Sources', f"{', '.join(run.query.sources)}"), field_line('Unique articles', f"{len(run.articles)}"), "", "SCIENTIFIC NOTE: Relevance scores are not evidence-quality scores.", ""]

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

            if a.journal:
                lines.append("JO  - " + a.journal)

            if a.publication_date:
                lines.append("PY  - " + a.publication_date[:4])

            if a.doi:
                lines.append("DO  - " + a.doi)

            if a.pmid:
                lines.append("AN  - PMID:" + a.pmid)

            if a.abstract:
                lines.append("AB  - " + a.abstract.replace("\n", " "))
            if a.metadata.get("date_validation") == "FUTURE_DATE_REVIEW":
                lines.append("N1  - REVIEW: FUTURE PUBLICATION DATE. " + a.metadata.get("date_validation_note", "Verify source date."))

            for url in a.source_urls:
                lines.append("UR  - " + url)

            lines.append("ER  - ")

            lines.append("")

        destination.write_text("\n".join(lines), encoding="utf-8")

    else:

        entries = []

        for index, a in enumerate(run.articles, 1):

            key = re.sub(r"[^A-Za-z0-9]", "", (a.authors[0].split()[-1] if a.authors else "article") + (a.publication_date or "")[:4] + str(index)) or f"article{index}"

            # LaTeX's special characters; an unescaped & or % in a title stops a
            # LaTeX run (Phase 7). Identifiers (DOI) keep all but the brace rules.
            latex_specials = {"&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
                              "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}

            def bib_escape(value: str, *, identifier: bool = False) -> str:
                escaped: list[str] = []
                for character in value:
                    if character == "\\":
                        escaped.append(r"\textbackslash{}")
                    elif character == "{":
                        escaped.append(r"\{")
                    elif character == "}":
                        escaped.append(r"\}")
                    elif character in latex_specials and not identifier:
                        escaped.append(latex_specials[character])
                    else:
                        escaped.append(character)
                return "".join(escaped)

            fields = [f"  title = {{{bib_escape(a.title)}}}"]

            if a.authors:
                fields.append("  author = {" + " and ".join(bib_escape(x) for x in a.authors) + "}")

            if a.journal:
                fields.append(f"  journal = {{{bib_escape(a.journal)}}}")

            if a.publication_date:
                fields.append(f"  year = {{{a.publication_date[:4]}}}")

            if a.doi:
                fields.append(f"  doi = {{{bib_escape(a.doi, identifier=True)}}}")

            if a.abstract:
                fields.append(f"  abstract = {{{bib_escape(a.abstract)}}}")
            if a.metadata.get("date_validation") == "FUTURE_DATE_REVIEW":
                fields.append("  note = {" + bib_escape("REVIEW: FUTURE PUBLICATION DATE. " + a.metadata.get("date_validation_note", "Verify source date.")) + "}")

            entries.append("@article{" + key + ",\n" + ",\n".join(fields) + "\n}")

        destination.write_text("\n\n".join(entries), encoding="utf-8")

    return destination
