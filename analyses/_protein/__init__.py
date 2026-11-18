"""Implementation of :mod:`analyses.protein` (split out in Phase 9).

Import from ``analyses.protein``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``constants``: alphabets, residue masses, hydropathy and flexibility
  scales, pKa values, the instability (DIWV) table
- ``sequence``: validation, resolved and positional views, sequence quality
- ``translation``: codon resolution (IUPAC), RNA translation
- ``composition``: amino-acid composition (``analyze_protein``), molecular weight, extinction
  coefficient, aromaticity, aliphatic and instability index
- ``charge``: net charge and isoelectric point
- ``profiles``: hydropathy, hydrophobic and transmembrane regions,
  secondary-structure propensity, flexibility, signal peptide
- ``report``: the combined analyses and the preview
"""
