"""Implementation of :mod:`analyses.mutation` (split out in Phase 8).

Import from ``analyses.mutation``; it re-exports every name defined here,
keeps ``__all__`` and is the only public interface. The modules, in
dependency order:

- ``sequence``: alphabets, IUPAC handling, normalization, translation
- ``models``: the result dataclasses and errors
- ``alignment``: global (banded) alignment, identity, edit distance, blocks
- ``detection``: mutation calls from an alignment, notation, left alignment
- ``coding``: codon and protein consequences, the coding-frame check
- ``hgvs``: HGVS c. and p. descriptions
- ``comparison``: comparing sequences and sets of variants, statistics
- ``context``: regions, sequence quality, context, spectrum, hotspots
- ``impact``: restriction-site, primer and guide overlap of mutations
- ``reporting``: validation, text and JSON reports, the advanced analysis
- ``selftest``: the built-in self-test
"""
