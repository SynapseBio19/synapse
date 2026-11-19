"""Implementation of :mod:`analyses.rna` (split out in Phase 9).

Import from ``analyses.rna``; it re-exports every name defined here and is
the only public interface. The modules, in dependency order:

- ``bases``: RNA base tables, codon table, transcription, counts, complements
- ``composition``: percentages, GC/AU content, ranges, skews, profiles,
  ambiguity, quality, molecular weight, complexity, homopolymers
- ``coding``: codon usage, Kozak and Shine-Dalgarno context, start/stop
  landscape, frame translation
- ``structure``: the hairpin screen (a heuristic; the energy model is
  ``analyses.rna_structure``)
- ``summary``: the feature score, the summary, ``analyze_rna``, the preview
- ``pipeline``: the RNA source of the session pipeline
"""
