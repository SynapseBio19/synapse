"""Implementation of :mod:`analyses.orf` (split out in Phase 9).

Import from ``analyses.orf``; it re-exports every name defined here and is
the only public interface. The modules, in dependency order:

- ``codes``: start/stop codons, policies and defaults, genetic-code tables,
  IUPAC codon handling, start-codon classification
- ``features``: per-ORF evaluation (ambiguity, codon usage, Kozak and
  Shine-Dalgarno context, CDS quality, heuristic score)
- ``finder``: ``find_orfs`` with partial and origin-spanning ORFs
- ``summary``: ``analyze_orfs``, density, frame statistics, relationships,
  ranking and comparison
- ``selection``: translation and analysis of a selected ORF,
  ``analyze_orfs_enhanced``
"""
