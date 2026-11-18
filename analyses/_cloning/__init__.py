"""Implementation of :mod:`analyses.cloning` (split out in Phase 8).

Import from ``analyses.cloning``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``cache``: the per-search result cache (``_strategy_search_cache``)
- ``sites``: sequences, enzyme data, overhangs, site finding, restriction maps
- ``ends``: cut ends, their compatibility and the junctions
- ``strategy``: evaluating and scoring one cloning strategy
- ``search``: searching all enzyme combinations, the automatic planner
- ``summary``: entry points, summaries, ranking and comparison
- ``selftest``: the built-in self-test

Tests that patch a helper must patch the module that calls it, for example
``analyses._cloning.sites._analyze_enzyme``; the cache state is
``analyses._cloning.cache._STRATEGY_SEARCH_CACHE``.
"""
