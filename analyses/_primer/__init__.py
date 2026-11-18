"""Implementation of :mod:`analyses.primer` (split out in Phase 7.1).

Import from ``analyses.primer``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``sequence``: DNA normalization, reverse complement, GC, repeats
- ``models``: the legacy design configuration and result dataclasses
- ``thermodynamics``: nearest-neighbor model, salt correction, Tm methods,
  condition presets
- ``structure``: hairpin, self-dimer and heterodimer screens
- ``specificity``: specificity screen against supplied references, binding sites
- ``design``: the legacy design engine (``design_primers``, ``analyze_primer``)
- ``professional``: the professional data model, thermodynamic scores,
  single-primer analysis and diagnostics
- ``screening``: the professional candidate pool and its Tm-tier screen
- ``ranking``: professional pair scoring, diversity selection, ranking
- ``workflow``: ``design_primers_professional`` and restriction-site design
- ``selftest``: the built-in self-tests and ``PRIMER_V3_PUBLIC_API``

The ``_v5``-``_v11`` prefixes are history, not versions: every such
function is live code (see the section comments).
"""
