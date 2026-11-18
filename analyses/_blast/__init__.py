"""Implementation of :mod:`analyses.blast` (split out in Phase 7.1).

Import from ``analyses.blast``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``paths``: the ``analyses`` directory and user-path resolution
- ``config``: endpoints, rate limits, programs, databases, local catalog
- ``models``: the dataclasses, ``BlastAnalysisError`` and hit ranking
- ``query``: query validation, FASTA loading, the ORF adapter
- ``parsing``: BLAST XML parsing
- ``ncbi``: the NCBI BLAST URL API client
- ``executables``: locating the BLAST+ programs
- ``runtime``: runtime alias views and the BLAST environment
- ``databases``: the local database manager
- ``local_engine``: local BLAST+ execution
- ``evidence``: evidence score, taxonomy grouping, alignment text
- ``reporting``: serialization, exports, text report, run comparison
- ``workflow``: ``run_blast_analysis``
- ``selftest``: the built-in self-test
"""
