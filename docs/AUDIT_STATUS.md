# Status of every audit recommendation (as of 2026-11-13)

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
| P1-1 | Test suite and golden master | IMPLEMENTED+TESTED | 109 tests at the end of the 2026-10 cycle; 314 at the end of Phase 6 (2026-11-02), including external reference tests against Biopython and primer3-py that skip without them, integration workflow tests, and interface and style enforcement tests; 13 golden snapshots; network and `input()` blocked in tests. |
| P1-2 | `SequenceRecord`, `Interval`, `Feature`, `AnalysisResult`, `store_result()`, no silent `except` | IMPLEMENTED+TESTED | Every module stores visibly; coordinate tests; failures are logged. |
| P1-3 | Mg²⁺ model, symmetry parameter, honest category labels | IMPLEMENTED+TESTED (external reference comparison) | von Ahsen 2001 implemented; known-answer test against the Primer3 manual example (ΔH, ΔS exact; Tm 50.550 vs 50.554 °C) and an independent re-implementation of the equations. The criterion "±1 °C against Biopython/Primer3 at three Mg levels" is met (Phase 6, 2026-11-02): against primer3-py 2.4.0 at Mg²⁺ 0, 1.5 and 3 mM, 2,916 comparisons agree within 0.0032 °C. This confirms the implementation against an independent one; it is not experimental validation of the model. |
| P1-4 | Restriction sites across the origin (circular) | IMPLEMENTED+TESTED | Implemented and unit tested (`wraps_origin`). Phase 6: the digest was compared with Biopython `Bio.Restriction` on pUC19 (132/132 agree); the map and the other site-based results did not honor topology (audit D3) and now do, also 132/132 against Biopython (2026-10-31). The NEBcutter comparison itself was not run; Biopython's REBASE-based engine served as the independent reference. |
| P1-5 | VCF: FILTER semantics, header cardinality A/R/G, PL/AD remapping, contig mapping, `.vcf.gz` | IMPLEMENTED+TESTED | FILTER semantics, header cardinality, PL/AD remapping and gzip implemented and tested. Contig-name matching added in Phase 6 (audit D12, 2026-10-31): exact by default, opt-in `chr_prefix` pairing of `1`/`chr1`, mitochondrial names never paired, and an unmatched contig is now flagged (`REFERENCE_CONTIG_NOT_FOUND`) instead of being skipped silently (`tests/test_phase6_tier1.py`). Not covered: contig-length checks and RefSeq accession names. |
| P1-6 | File selection (FASTA/GenBank), topology, SHA-256 | IMPLEMENTED+TESTED | Own files load; hashes recorded in every result. |
| P1-7 | BLAST: effective parameters, version, DB date, privacy notice, RID resume | PARTIAL | Reported parameters, database provenance and the privacy confirmation are implemented. RID resume is implemented (Phase 6, P1-7): `NCBIBlastClient.resume()`, a menu entry and `core.api.blast_resume()`, tested against a stubbed client rather than the live service. Phase 6: local BLAST+ runs are tested against the bundled `data/local_test` database and small `makeblastdb` databases (blastn, blastp), including the programmatic database-path entry (D6, S7). |
| P1-8 | Literature: date formats, throttling, backoff, optional API key | IMPLEMENTED+TESTED | Fixtures for all date formats; retry test with a simulated HTTP 503 and `Retry-After`. |
| P1-9 | ORF: genetic codes, start policy, organism mode, 5′-partial, circular | COMPLETE (not externally validated) | Genetic codes 1, 2, 4, 11 and the start policies; the initiator set follows the genetic code (D1, D2); the RNA codon view uses the ORF's code (2026-10-31). Since 2026-11-11 (`f020a1d`): ORFs across the origin of circular sequences and optional 5′-partial ORFs, with a rotation-invariance test. Organism mode is covered by genetic code plus start policy, without a separate switch. Open validation: no ORFfinder comparison (no native Windows build; needs WSL or Linux). |
| P1-10 | Mutation: CDS start/frame from ORF, frame warning, length warning | IMPLEMENTED+TESTED | CDS from the selected ORF (both strands), plausibility warnings, 2 kb note, "HGVS-like" labeling. |
| P1-11 | `pyproject.toml`, Python version, Git, logging, README with methods | IMPLEMENTED+TESTED | Git history from a Stage 0 baseline; rotating log file; `docs/METHODS.md`. |
| P1-12 | Reproducible analysis report | IMPLEMENTED+TESTED | Menu `[R]`; JSON (complete) and Markdown; input hash, parameters, versions, limitations; end-to-end menu test. |

## P2 – structure, breadth, performance

| ID | Subject | Status | Evidence / open points |
|---|---|---|---|
| P2-1 | Consolidate the primer patch chain behavior-preserving | PARTIAL (by decision) | 13 unreachable superseded definitions removed in the previous cycle; in Phase 6 (D16) only the last binding of each name was found reachable and the 24 private definitions kept alive by superseded aliases were removed (129 -> 105 module-level callables, every survivor bytecode-identical). Since then each role has exactly one live implementation; the remaining `_v5`–`_v11` prefixes are historical names, not a wrapper chain. Renaming them would change about 40 references without changing behavior, so it was not done (decision of 2026-11-12). The primer candidate screen was rebuilt for speed with proven identical results (2026-11-12, `aba3fb0`). |
| P2-2 | Separate UI from compute modules; `core/seqtools.py` against duplication | COMPLETE | `core/api.py` offers every analysis without `input()`; genetic-code and IUPAC tables are shared via `core/seqtools.py`. Since 2026-11-12 (`50dea66`) the dialogs of all eleven analysis modules are in `ui/<module>_menu.py` and no analysis module imports `ui` (outside `if __name__ == "__main__"` script blocks); string layout comes from `core/textformat.py`. Evidence: 55 rendered screens byte-identical before and after; `LayeringTests` runs every `core.api` analysis with `ui` unimportable. |
| P2-3 | CRISPR: cut position, U6 rules, local CFD scoring, exclude the on-target site | DEFERRED / partly REJECTED | Phase 6: the modeled cut site is reported with its source (audit D8) and the on-target site is excluded from the local specificity screen (audit D9, 2026-10-31). The rest was not implemented in this session. **CFD scoring is rejected** for now: it requires the published mismatch weight matrix, and reconstructing those numbers from memory would mean inventing scientific data. It can be added once the matrix is supplied as a data file with a citation. |
| P2-4 | Cross-checks (primer/guide against variants, REF vs. sequence) | COMPLETE | REF vs. reference sequence: `check_reference_alleles` (P0-7). Primer and guide binding sites against variants: `core.api.cross_check_variants` (2026-11-11, `2c68e61`; 3′ end of primers, seed and PAM of guides as critical parts) and, since 2026-11-13 (`cdbf69c`), in the primer and CRISPR menus for the variants of a VCF contig mapped onto the loaded sequence. |
| P2-5 | Honest labels; preserve "?" positions; G·U / RNAfold; rename complexity | PARTIAL | Unresolved residues keep their position (Stage 5); labels and validation status are documented in `docs/METHODS.md`; heuristic scores are no longer called confidence (D18); "complexity" is documented as normalized Shannon entropy that does not detect repeats (2026-11-02). G·U pairs are counted separately in the RNA hairpin screen since 2026-11-11 (`fedb750`). Open: RNAfold (an external binary, not installed) and the rename of the `complexity` field (kept for compatibility). |
| P2-6 | Performance (banded alignment, sliding CpG count, primer pre-filter) | COMPLETE for the 10 kb target | Recognition-site scan cache (S18); cloning planner 21.7 s -> 3.2 s (`b4b16bb`); mutation alignment one byte per cell and an exact band for similar sequences (10 kb 23 s -> 0.18 s); banded Levenshtein similarity (10 kb 0.10 s); primer design 10 kb 135.6 s -> 25.6 s (`15e2ff1`) -> about 10 s CPU time under heavy background load, 71 s for the previous code under the same load (`aba3fb0`). Every change has identical results by construction and by comparison. Remaining: unrelated long sequences still align in O(n·m) time (about 23 s at 10 kb). |
| P2-7 | Type IIS enzymes on both strands, methylation flags | COMPLETE (definitions and digests) | Methylation notes (ClaI, XbaI, DpnI, HpaII). Since 2026-11-11 (`4939688`): BsaI, BsmBI, BbsI, SapI on both strands, cutting outside their sites; definitions 4/4 and 140 cuts vs Biopython 1.88, 0 differences. Not part of the cloning planner; Golden Gate planning not implemented. |

## P3 – long-term

| ID | Subject | Status |
|---|---|---|
| P3-1 | Evidence-bound experiment planner | DEFERRED (depends on P1-2, P1-12, both now in place). The menu entry says "not implemented" (Phase 6, D15). |
| P3-2 | Primer/guide specificity against genomes (BLAST, Cas-OFFinder) | DEFERRED (needs local databases) |
| P3-3 | Web or desktop front end | DEFERRED (depends on P2-2; `core/api.py` is the entry layer for it) |

## Validation that remains open

These were not possible in the development environment and are the
recommended next steps:

1. `bcftools norm` comparison for variant normalization (P0-7).
2. ~~Primer3 / Biopython Tm comparison at several Mg²⁺ levels (P1-3).~~
   Done 2026-11-02 with primer3-py 2.4.0: 2,916 comparisons over Mg²⁺
   0-3 mM, max |ΔTm| 0.0032 °C (`docs/METHODS.md`).
3. NEBcutter comparison for a circular plasmid digest (P1-4).
4. ORFfinder comparison on reference sequences (P1-9).
5. ~~A systematic REBASE re-check of all enzyme entries (P0-4).~~ Done
   2026-11-13 against REBASE version 611 (`emboss_e.611`): 26 of 26
   definitions identical (`docs/METHODS.md`).
6. ~~Local BLAST+ runs~~ done in Phase 6 (BLAST+ 2.17.0 on the bundled
   database and on `makeblastdb` databases, blastn and blastp).
   ~~One live remote NCBI run~~ done 2026-11-13 (pUC19 L09137.2 bases
   401-700 against `nt`; `docs/METHODS.md`); the first attempt found D25.
7. ~~An independent review of the primer ranking chain~~ done in Phase 6
   (D16 part 2, reachability review; P2-1). Renaming the surviving
   `_v5`–`_v11` functions was decided against (2026-11-12, P2-1 above).
8. ~~A Python 3.10 suite run~~ done 2026-11-12 with the official embeddable
   CPython 3.10.11 in a project-local folder: 452 tests, all passing
   (`docs/TOOL_AND_DATA_LEDGER.md`).
9. The instability-index weights against the Guruprasad et al. (1990)
   table itself; Phase 6 corrected the IG weight against Biopython's
   table (`ProtParamData.DIWV`, which cites Guruprasad et al. 1990),
   a single secondary source.

## Phase 6 audit findings

The Phase 6 initial audit (`docs/PHASE6_INITIAL_AUDIT.md`) recorded 19
confirmed defects (D1-D19) and 14 suspected or unverified concerns, a
list that grew to 24 (S1-S24) as the work found more. Their status at the
end of Phase 6:

- **D1-D19**: all corrected; D18 (the "confidence" vocabulary) in Stage G
  with a documented golden-master change.
- **S1-S24**: all assessed. Confirmed ones are corrected; S1, S2 are not
  defects (verified by execution); S24 is measured with its residual risk
  recorded; **S11** (live literature providers) was checked live on
  2026-11-12 under the approval policy adopted that day (`docs/METHODS.md`).
- One further defect was found in Stage I by the new reference tests: the
  instability weight for IG (corrected; `docs/CHANGES.md`).
- The completion run of 2026-11-12/13 found and corrected six more, each
  with a test that fails on the previous code (`docs/CHANGES.md`): D20 the
  primer structural penalty was always zero; D21 `tm_method`
  "nearest-neighbor" crashed candidate generation; D22 an ORF without a
  sequence crashed `orf_reference_cds`; D23 the primer pair risk assessment
  never reached exports; D24 the BLAST text export and query file were
  written with CR CR LF; D25 a brief network error while polling lost a
  submitted NCBI search and its RID. Two BLAST dialog faults (hit
  selector, query-type question) were corrected as Stage F items.

Each item's evidence, commit and validation level is in
`docs/PHASE6_PROGRESS.md`; the per-feature status is in
`docs/FEATURE_STATUS.md`.
