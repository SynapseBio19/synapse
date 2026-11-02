# Status of every audit recommendation (as of 2026-10-30)

Assessment categories:

- **IMPLEMENTED+TESTED** – implemented and covered by a known-answer or
  regression test in `tests/`; where the audit named an acceptance
  criterion, that criterion is checked in the test suite.
- **IMPLEMENTED (not independently validated)** – implemented and unit
  tested, but the audit's acceptance criterion required a comparison with an
  external tool or data set that was not available in this environment.
- **PARTIAL** – the substance is implemented, a named sub-item is not.
- **DEFERRED** – not started; reason given.
- **REJECTED** – deliberately not implemented; reason given.

No item below is claimed to be experimentally validated. Synapse has not
been benchmarked against established tools except where a specific check is
named.

## P0 – scientific defects

| ID | Subject | Status | Evidence / open points |
|---|---|---|---|
| P0-1 | BLAST coverage = union of HSP intervals; evidence band only at E ≤ 1e‑3 | IMPLEMENTED+TESTED | Acceptance criterion met: the values of the audited run (32 bp HSP on a 99 bp query) recompute to 32.32 % coverage, and E = 7.45 gives `NOT_SIGNIFICANT` (`test_p0_regressions.py`). The stored file of that run, `blast_results/TEST_8_REFERENCE_blast.json`, predates the fix and is kept only as a labeled legacy record (Phase 6, D7). Threshold configurable. |
| P0-2 | Shared ORF adapter for CRISPR and BLAST (0→1-based, strand) | IMPLEMENTED+TESTED | `core/orf_context.py`; tests check that every reported guide position matches the original DNA on both strands. |
| P0-3 | End compatibility checks overhang polarity; overhang derived from cut coordinates | IMPLEMENTED+TESTED | Fallback table removed, `get_overhang` derived from cut positions; pairwise tests incl. identical sequence with opposite polarity. The audit's criterion "pair matrix equals the NEB compatibility table" was checked only for the enzymes in the database, not against a copy of the NEB table. |
| P0-4 | Enzyme data (ClaI, NdeI) with source and date | PARTIAL | ClaI corrected (AT^CGAT, 2 nt 5′ CG), NdeI overhang TA confirmed, `ENZYME_DATABASE_VERSION` added. A systematic re-check of all 22 entries against REBASE was not possible (no access from this environment); the entries carry a reviewed-date comment instead. |
| P0-5 | Ambiguity-dependent (POSSIBLE) sites excluded from digests | IMPLEMENTED+TESTED | `include_possible=False` by default, reported separately; test: 110 bp with 10 N → 0 cuts, 1 fragment. Phase 6 (audit D4, 2026-10-31): the cloning module had its own matcher without this distinction; it now uses the restriction matcher, counts only CONFIRMED sites, and does not treat an enzyme with POSSIBLE sites as a unique cutter (`tests/test_phase6_tier1.py::CloningAmbiguousSiteTests`). |
| P0-6 | Start-loss only for codon 1 | IMPLEMENTED+TESTED | Internal ATG→ATA = 1 missense, 0 start-loss. |
| P0-7 | Variant normalization per Tan et al. 2015, REF check | IMPLEMENTED (not independently validated) | Algorithm, REF check and `REF_MISMATCH` flag implemented and unit tested. The acceptance criterion "result equals `bcftools norm`" could not be run (tool not installable here). Recommended as the next validation step. |
| P0-8 | Initiator methionine at an ORF's first codon | IMPLEMENTED+TESTED | GTG-ORF → "MKF"; configurable, bound to the genetic code. Phase 6 (audit D1, 2026-10-31): the initiator set now comes from the genetic code (NCBI starts, plus table-1 GTG policy), so CTG (table 1) is also decoded as Met, and the protein/API path no longer turns a non-initiation first codon into Met; checked against Biopython for tables 1, 2, 4, 11. |
| P0-9 | ORF defaults: ATG only, 75 nt, longest per stop, length ranking | IMPLEMENTED+TESTED | Approved by the user before the change. Legacy settings selectable and proven equivalent to the old behavior by a golden-master test. |

## P1 – robustness, reproducibility, expected options

| ID | Subject | Status | Evidence / open points |
|---|---|---|---|
| P1-1 | Test suite and golden master | IMPLEMENTED+TESTED | 109 tests (`unittest`, standard library); 12 golden scenarios; network and `input()` blocked in tests. |
| P1-2 | `SequenceRecord`, `Interval`, `Feature`, `AnalysisResult`, `store_result()`, no silent `except` | IMPLEMENTED+TESTED | Every module stores visibly; coordinate tests; failures are logged. |
| P1-3 | Mg²⁺ model, symmetry parameter, honest category labels | IMPLEMENTED (not independently validated) | von Ahsen 2001 implemented; known-answer test against the Primer3 manual example (ΔH, ΔS exact; Tm 50.550 vs 50.554 °C) and an independent re-implementation of the equations. The criterion "±1 °C against Biopython/Primer3 at three Mg levels" was not run: neither package could be installed or executed here. |
| P1-4 | Restriction sites across the origin (circular) | IMPLEMENTED+TESTED | Implemented and unit tested (`wraps_origin`). Phase 6: the digest was compared with Biopython `Bio.Restriction` on pUC19 (132/132 agree); the map and the other site-based results did not honor topology (audit D3) and now do, also 132/132 against Biopython (2026-10-31). The NEBcutter comparison itself was not run; Biopython's REBASE-based engine served as the independent reference. |
| P1-5 | VCF: FILTER semantics, header cardinality A/R/G, PL/AD remapping, contig mapping, `.vcf.gz` | IMPLEMENTED+TESTED | FILTER semantics, header cardinality, PL/AD remapping and gzip implemented and tested. Contig-name matching added in Phase 6 (audit D12, 2026-10-31): exact by default, opt-in `chr_prefix` pairing of `1`/`chr1`, mitochondrial names never paired, and an unmatched contig is now flagged (`REFERENCE_CONTIG_NOT_FOUND`) instead of being skipped silently (`tests/test_phase6_tier1.py`). Not covered: contig-length checks and RefSeq accession names. |
| P1-6 | File selection (FASTA/GenBank), topology, SHA-256 | IMPLEMENTED+TESTED | Own files load; hashes recorded in every result. |
| P1-7 | BLAST: effective parameters, version, DB date, privacy notice, RID resume | PARTIAL | Reported parameters, database provenance and the privacy confirmation are implemented. RID resume is implemented (Phase 6, P1-7): `NCBIBlastClient.resume()`, a menu entry and `core.api.blast_resume()`, tested against a stubbed client rather than the live service. Phase 6: local BLAST+ runs are tested against the bundled `data/local_test` database and small `makeblastdb` databases (blastn, blastp), including the programmatic database-path entry (D6, S7). |
| P1-8 | Literature: date formats, throttling, backoff, optional API key | IMPLEMENTED+TESTED | Fixtures for all date formats; retry test with a simulated HTTP 503 and `Retry-After`. |
| P1-9 | ORF: genetic codes, start policy, organism mode, 5′-partial, circular | PARTIAL | Genetic codes (1, 2, 4, 11) and the start policies are implemented and tested; since Phase 6 (audits D1, D2, 2026-10-31) the initiator set follows the genetic code and an ORF is translated with its own table in the ORF and protein/API paths (and since Phase 6 the RNA module's codon view uses the selected ORF's table too). **Organism mode, 5′-partial ORFs and ORFs across the origin are deferred.** The criterion "comparison with ORFfinder on reference sequences" was not run (no access). |
| P1-10 | Mutation: CDS start/frame from ORF, frame warning, length warning | IMPLEMENTED+TESTED | CDS from the selected ORF (both strands), plausibility warnings, 2 kb note, "HGVS-like" labeling. |
| P1-11 | `pyproject.toml`, Python version, Git, logging, README with methods | IMPLEMENTED+TESTED | Git history from a Stage 0 baseline; rotating log file; `docs/METHODS.md`. |
| P1-12 | Reproducible analysis report | IMPLEMENTED+TESTED | Menu `[R]`; JSON (complete) and Markdown; input hash, parameters, versions, limitations; end-to-end menu test. |

## P2 – structure, breadth, performance

| ID | Subject | Status | Evidence / open points |
|---|---|---|---|
| P2-1 | Consolidate the primer patch chain behavior-preserving | PARTIAL | 13 unreachable superseded definitions (1,014 lines) removed; equivalence shown by comparing the compiled bytecode of all live functions and classes before and after, plus the golden master. The live v5–v11 wrapper chain for ranking was **kept**: flattening it would change executed code and needs its own stage. Phase 6 (audit D16, 2026-10-31): the last duplicated module-level name was removed by defining the Tm model under its private name instead of shadowing the public one; equivalence shown again by bytecode fingerprinting of all 129 module-level callables. Phase 6 part 2: the reachability review recorded as open validation item 7 was carried out (AST binding order, project-wide reference scan, `sys.settrace` on a real design run). It showed that only the last binding of each public name was ever reachable, so the 7 superseded aliases and the 24 private definitions they alone kept alive were removed; 129 -> 105 module-level callables with every survivor bytecode-identical. Merging the surviving generations is deliberately not done. |
| P2-2 | Separate UI from compute modules; `core/seqtools.py` against duplication | PARTIAL | `core/api.py` offers every analysis without `input()`; genetic-code and IUPAC tables are shared via `core/seqtools.py`. The existing menus still live inside the analysis modules, so the criterion "core modules without `input()`" is met by the new API layer but not by the modules themselves. |
| P2-3 | CRISPR: cut position, U6 rules, local CFD scoring, exclude the on-target site | DEFERRED / partly REJECTED | Phase 6: the modeled cut site is reported with its source (audit D8) and the on-target site is excluded from the local specificity screen (audit D9, 2026-10-31). The rest was not implemented in this session. **CFD scoring is rejected** for now: it requires the published mismatch weight matrix, and reconstructing those numbers from memory would mean inventing scientific data. It can be added once the matrix is supplied as a data file with a citation. |
| P2-4 | Cross-checks (primer/guide against variants, REF vs. sequence) | PARTIAL | REF vs. reference sequence is implemented (`check_reference_alleles`, P0-7). Primer/guide-versus-variant consistency checks are deferred. |
| P2-5 | Honest labels; preserve "?" positions; G·U / RNAfold; rename complexity | PARTIAL | Unresolved residues keep their position in all position-dependent protein analyses (Stage 5, tested); labels and validation status are documented in `docs/METHODS.md`; the ORF "Confidence" label is now "Heuristic score", and since 2026-11-02 (D18) the data model and exports say `heuristic_score` as well. G·U pairing, optional RNAfold and the complexity renaming are deferred. |
| P2-6 | Performance (banded alignment, sliding CpG count, primer pre-filter) | DEFERRED | Only a length warning above 2 kb was added (P1-10). The criterion "10 kb in < 10 s per module" is not met for the mutation module. |
| P2-7 | Type IIS enzymes on both strands, methylation flags | PARTIAL | Methylation notes are implemented (ClaI, XbaI, DpnI, HpaII). **Type IIS support is deferred**; the forward-strand-only search remains correct for the palindromic enzymes in the database, and no Type IIS enzyme was added, so no wrong result can arise from this today. |

## P3 – long-term

| ID | Subject | Status |
|---|---|---|
| P3-1 | Evidence-bound experiment planner | DEFERRED (depends on P1-2, P1-12, both now in place) |
| P3-2 | Primer/guide specificity against genomes (BLAST, Cas-OFFinder) | DEFERRED (needs local databases) |
| P3-3 | Web or desktop front end | DEFERRED (depends on P2-2; `core/api.py` is the entry layer for it) |

## Validation that remains open

These were not possible in the development environment and are the
recommended next steps:

1. `bcftools norm` comparison for variant normalization (P0-7).
2. Primer3 / Biopython Tm comparison at several Mg²⁺ levels (P1-3).
3. NEBcutter comparison for a circular plasmid digest (P1-4).
4. ORFfinder comparison on reference sequences (P1-9).
5. A systematic REBASE re-check of all enzyme entries (P0-4).
6. Local BLAST+ runs against `data/blast_db/` and one live remote run (P1-7).
7. An independent review of the primer ranking chain before it is flattened
   (P2-1).
