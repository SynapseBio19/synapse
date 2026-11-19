"""Implementation of :mod:`analyses.sequence` (split out in Phase 9).

Import from ``analyses.sequence``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``bases``: base tables, normalization, counts, complement
- ``composition``: percentages, GC/AT content and ranges, skews, quality,
  molecular weight
- ``profiles``: sliding-window GC and complexity profiles, homopolymers
- ``cpg``: CpG statistics and CpG islands
- ``summary``: sequence forms, interpretation, the summary and
  ``analyze_sequence``
- ``selftest``: the built-in self-test
"""
