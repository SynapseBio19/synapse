"""Implementation of :mod:`analyses.restriction` (split out in Phase 9).

Import from ``analyses.restriction``; it re-exports every name defined here
and is the only public interface. The modules, in dependency order:

- ``sites``: sequence validation, IUPAC site matching (also Type IIS),
  cut positions, overhangs, per-site and per-enzyme analysis
- ``digest``: restriction maps, cut positions, linear and circular
  fragments, digests, gel bands, cutter lists
- ``compatibility``: sticky ends, cloning compatibility, enzyme-pair
  recommendations, double digests
- ``api``: the combined analysis, summaries and the ``get_*`` / ``run_*``
  entry points, the enzyme-database check
- ``selftest``: the built-in self-test
"""
