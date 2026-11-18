"""Implementation of :mod:`analyses.variant_analysis` (split out in Phase 9).

Import from ``analyses.variant_analysis``; it re-exports every name defined
here, keeps ``__all__`` and is the only public interface. The modules, in
dependency order:

- ``models``: constants, the dataclasses, ``VariantAnalysisError`` and
  ``classify_variant_type`` (used by ``Variant``)
- ``vcf``: VCF parsing (header numbers, INFO, genotypes, ``VariantRecord``
  and the multi-allelic split)
- ``normalization``: trimming, left normalization, contig matching,
  reference checks
- ``annotation``: annotation from INFO and from the mutation analysis,
  evidence, prioritization, feature overlap
- ``statistics``: filters, statistics, cohort summary, set comparison,
  grouping and intervals
- ``export``: VCF and TSV export, the text report
- ``workflow``: ``analyze_variants`` and ``analyze_vcf_text``
- ``selftest``: the built-in self-test
"""
