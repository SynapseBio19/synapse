"""Implementation of :mod:`analyses.literature` (split out in Phase 10).

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
