# Phase 6 progress record

Running record of the Phase 6 programme. The audit itself is in
`docs/PHASE6_INITIAL_AUDIT.md`; the feature matrix in
`docs/FEATURE_STATUS.md`. This file says what has been done, how it was
checked, and what is still open. It is updated with every work package.

Validation vocabulary used below:

- **internal test** – a unit, regression or golden-master test in `tests/`.
- **external reference comparison** – a test that compares Synapse with an
  independent implementation (Biopython, primer3-py); skipped when the
  package is absent.
- **external tool execution** – a run of a third-party program (BLAST+,
  bcftools) on real input.
- **literature support** – a cited published convention or data source.

## Stage A – revalidation (2026-10-31)

State reconstructed from the repository, not from earlier summaries.

| Item | Verified state |
|---|---|
| Branch / tree | `main`, clean at `ab4bd2a` |
| Phase 6 commits | `5242fcb`, `5d8c73a`, `824ca5f` (housekeeping, audit); `d4544ad` D12; `5c098e3` D4; `47146b0fcead81a43906881abf2b02bc2a9991ec` D3; `7334296` D5; `dc780f6` D1; `ab4bd2a` D2 |
| Hash discrepancy | The D3 hash quoted in the revalidation prompt (`…abf2f02bc…`) does not exist; the actual commit is `…abf2b02bc…` (same subject, same short hash `47146b0`). |
| Tests | 151/151 pass (`python -m unittest discover -s tests -t .`) |
| pyright | 1.1.414, `--pythonversion 3.10`, 49 files, 134 errors, 0 warnings (pre-existing; 49 in `analyses/primer.py`). Every work package is checked for *new* diagnostics against this set. |
| Python | venv 3.14.6; target compatibility 3.10 (pyright) |
| Optional packages | biopython 1.88, primer3-py 2.4.0 installed in the venv (test-only) |
| External tools | BLAST+ 2.17.0 (`blastn`, `makeblastdb`, `blastdbcmd`) present; `bcftools` not present |
| Skills | 11 local skills under `.claude/skills/` (ignored by Git) |
| MCP | PubMed, bioRxiv, Context7 available for literature and documentation checks |

Confirmed again by a probe of every module with DNA, RNA, protein,
IUPAC, gapped and invalid input (D14):

- BLAST classifies RNA and a T/U mixture as protein and then rejects `U`
  as an invalid protein character; a protein made only of letters that are
  also IUPAC nucleotide codes (e.g. `MKVHDANWYSK`) is classified as DNA.
- BLAST's program table uses `"protein"`, its queries `"PROTEIN"`: local
  `blastp`/`tblastn` always fail the query-type check, and the menu offers
  no program at all for a protein query (D14-F2).
- `SequenceRecord.alphabet` is always `"DNA"`; neither loader reads a
  molecule type (D14-F3).
- The main input gate accepts the empty string as valid DNA.

## Work packages

| # | Finding | Commit | Validation | Status |
|---|---|---|---|---|
| 1 | D14-F2 BLAST query-type vocabulary; protein BLAST usable | `b8f5ea8` | internal tests; external tool execution (BLAST+ 2.17.0 `makeblastdb` + `blastp`) | Corrected; suite 156/156; pyright unchanged |
| 2 | Central alphabets and classifier (F12/D13 groundwork) | `d057ac1` | internal tests; external reference comparison (`Bio.Data.IUPACData`); literature (PMIDs 2582368, 6692818) | Added, not yet wired in; suite 167/167; pyright unchanged |
| 3 | D14-F3, F11 molecule type from GenBank/FASTA; type in `SequenceRecord`, provenance and report | `776de4a` | internal tests; NCBI `[moltype=]` modifier checked in NCBI's FASTA defline documentation | Corrected; suite 173/173; pyright unchanged |
| 4 | D14-F1 BLAST query type from the classifier; S20 J/U/O accepted | `7c716a2` | internal tests; external tool execution (blastp on the IUPAC-letter peptide; BLAST+ acceptance of U, O, J checked directly) | Corrected; suite 180/180; pyright unchanged |
| 5 | D13 + D14-F4/F5/F7: one input contract; menu type gate with user declaration; API type checks; protein record route | `162ca16` | internal tests; end-to-end `main.py` run with protein, undetermined, RNA and invalid records | Corrected; suite 189/189; pyright unchanged |
| 6 | D14-F6 protein codes B/Z/J/X/U/O keep position; masses and ranges | `343ffb7` | internal tests; external reference comparison (Biopython `protein_weights`, `molecular_weight`); one Stage 5 assertion updated (encoded the deletion) | Corrected; suite 194/194; pyright unchanged |
| 7 | D14-F8/F9/F10 ORF search, transcription and GC counting check/normalize input | `d40d328` | internal tests; probe re-run (no silent results for foreign types) | Corrected; D14 and D13 complete; suite 197/197; pyright unchanged |
| 8 | D6 + S7 programmatic local BLAST by database prefix, with provenance | `5b120ec` | internal tests; external tool execution (bundled database, blastn + blastdbcmd) | Corrected; suite 200/200; pyright unchanged |
| 9 | D7 stale BLAST artifact labeled as legacy; AUDIT_STATUS cites the test | `e1a48c6` | documentation; the README's numbers match `tests/test_p0_regressions.py` | Corrected (file kept; deletion is the owner's decision) |
| 10 | D9 CRISPR on-target site excluded from the local specificity screen | `3514617` | internal tests (both strands, identical copy, near match) | Corrected; suite 203/203; pyright unchanged; S23 recorded |
| 11 | D8 CRISPR cut site per system, with citations | `ca2ec39` | internal tests (both strands, all four systems); literature (PMIDs 22745249, 25830891, 30166441, 26422227) | Corrected; suite 206/206; pyright unchanged |
| 12 | D10 dead CRISPR config fields act; D11 exports in `exports/crispr` | `8145556` | internal tests (defaults unchanged) | Corrected; suite 210/210; pyright unchanged |
| 13 | D19 pipeline `run_orf` search settings; fixed `completed` shape | `e47cb85` | internal tests (equality with `core.api.find_orfs`) | Corrected; suite 212/212; pyright unchanged |
| 14 | D15 planner entry labeled "not implemented" | `f3a09c1` | end-to-end `main.py` test | Corrected (planner itself remains P3-1, deferred); suite 213/213 |
| 15 | S15 variant-to-mutation bridge works; S16 REF_MISMATCH no longer leaks into the input | `c71963b` | internal tests; pyright 134 → 132 (the two S15 diagnostics) | Corrected; suite 215/215 |
| 16 | S17 cloning digest uses the top-strand cut | `c743beb` | internal tests; agreement with `restriction.perform_digest` for all 231 enzyme pairs, linear and circular | Corrected; suite 217/217 |
| 17 | S19 start-codon classification and scores use the genetic code | `6c20e69` | internal tests; golden master unchanged | Corrected; suite 219/219 |
| 18 | S23 PAM-proximal guide seed; S3 real RTT coverage check; S4/S5 thresholds documented and sourced | this package | internal tests; literature verified in PubMed (PMIDs 24476820, 27096365, 29160308, 31634902, 33811918) | Corrected; suite 228/228; pyright unchanged. Done by Opus 5 on a restricted assignment after a safeguard interrupted the Opus 5.5 session. |
| 19 | D16 part 1: one definition per name in `analyses/primer.py` | this package | bytecode fingerprint of all 129 module-level callables unchanged; pyright `reportRedeclaration` gone | Partially corrected; v5–v11 chain deferred; suite 229/229 |
| 20 | S21, S10 GenBank reader reports unread symbols and locations; S8 failing `to_dict()` logged; S22 variant reference must be DNA | this package | internal tests; bundled GenBank fixture unaffected | Corrected; suite 236/236 |
| 21 | S1, S2 assessed as not defects (locked by tests); S9 quantified by measurement; S11, S12 resolved in documentation | this package | execution on constructed minus-strand cases; `tracemalloc` measurement at three sizes | Assessed; suite 240/240 |
| 22 | D16 part 2: superseded primer generations removed after a reachability review | this package | AST binding order + project-wide reference scan + `sys.settrace` on a real design run; bytecode identical for all 105 surviving callables; pyright 134 -> 112 | Corrected; 1,326 lines removed; suite 240/240 |

| 23 | P1-9 RNA codon view uses the selected ORF's genetic code | this package | internal tests; +1 frame preview equals `translate_orf_sequence` for a table-2 ORF | Corrected; suite 245/245 |

| 24 | S18 recognition-site scans cached; planner no longer rescans | this package | cProfile before (5,002 scans for 16 combinations); measured 11.8 s -> 0.24 s at 8 enzymes, 100-110 s -> 0.57 s at 10 | Corrected; suite 248/248 |

| 25 | S14 API covers cloning, variant filtering, mutation reports, BLAST and literature | this package | internal tests incl. both privacy gates; local BLAST+ run | Resolved for the API; menus stay with Stage F (P2-2); suite 256/256 |

| 26 | S6/P1-7 resume an NCBI search from its request ID | this package | internal tests against a stubbed client (asserts no submission); no live NCBI call | Implemented; suite 265/265 |

| 27 | S13 session moved into `main()` behind a `__main__` guard | this package | `git diff -w` shows only the 12 added lines; fresh-interpreter import test | Corrected; suite 267/267; S24 recorded |

| 28 | Stage F: one separator width through `ui/style.py` (11 -> 1); `[0]` means Back/Cancel/Exit only (7 -> 3); uniform message prefixes; `TERMINAL_UI_STANDARD.md`; enforcement tests | `5db21ab` + this package | interface inventory before and after; enforcement tests verified against deliberate violations | D17 corrected; suite 277/277 |
| 29 | Stage F: BLAST menu renumbered so SEARCH reads 1-4 instead of 1, 2, 3, 14; numbering rule documented and enforced | this package | each branch checked against its new label; contiguity test for the BLAST and main menus | Corrected; suite 279/279 |
| 30 | D18 heuristic scores renamed from "confidence" in data, exports and messages (ORF, signal-peptide screen, cloning score component, primer warning); golden key renamed with values proven identical | this package | golden before/after hashes; key-mapped-back equality; unmodified golden test; new tests fail on the previous code | Corrected; suite 288/288 |
| 31 | D18 part 2: the sequence-type classifier reports an `evidence` level (CONFIRMED/PROBABLE/AMBIGUOUS/INVALID), not a "confidence"; source check extended to `core/`, `fileio/`, `ui/`, `main.py` | this package | the extended checks fail on the previous code; values and rules unchanged; no golden snapshot affected | Corrected; suite 289/289 |
| 32 | Stage G: one notation per unit in printed output (%, °C, cal/(mol·K), Da); primer positions and target input 1-based inclusive like the other menus | this package | inventory of every f-string literal; source check fails on the previous code; primer coordinates verified by slicing | Corrected; suite 292/292 |
| 33 | Stage G (6.4): payload warnings (primer, CRISPR) reach the exported summary and the Markdown report on every path | this package | measurement over every offline `core.api` analysis; new tests fail on the previous code | Corrected; suite 294/294; menu-path parameters for sequence/restriction/cloning recorded as open |
| 34 | Stage G (6.4): menu-path restriction results record enzymes and topology; restriction and cloning results record the enzyme data version on both paths | this package | new test fails on the previous code | Corrected; suite 295/295; menu-path sequence-overview and cloning input parameters still open |
| 35 | Stage H editorial review part 1: American spelling throughout (92 words in docs, 6 code text sites); one wording for recurring messages; tone and precision fixes; enforcement test | this package | inventory of every string literal and document; spelling test fails on the previous text | Done; suite 296/296 |
| 36 | BLAST E-values with three significant figures everywhere; `n/a` when absent | `4e5dbdc` | helper assertions (the sample XML cannot distinguish) | Corrected |
| 37 | `docs/STYLE_GUIDE.md` (9.1-9.6) with an enforcement table checked by the suite | this package | guide-consistency test verified with a deliberate violation | Done |

## Paused (2026-10-31)

Planned pause at the owner's request; see `docs/PHASE6_HANDOVER.md` for the
exact state, the open findings and the continuation steps. Stages F–I have
not started.

## Open after the current package

D17 and D18 (stages F and G), plus S6/P1-7, S13, S14/P2-2, S18 and P1-9.
Corrected or assessed: D1–D16, D19, S1–S5, S7–S12, S15–S17 and S19–S23.
S11 stays open pending the owner's approval for live external requests.

## Classification of the outstanding audit items (2026-11-01)

Every item is classified against the original Master Prompt, with the
reason it stands where it does. None is silently closed.

| Item | Master Prompt reference | State | Why |
|---|---|---|---|
| **P0-4** enzyme definitions vs REBASE | 3.6 "Verify enzyme definitions against authoritative sources" | **Partially done; blocked for the rest** | ClaI was corrected and NdeI confirmed in the previous cycle, and `ENZYME_DATABASE_VERSION` records the review date. A systematic re-check of all 22 entries needs REBASE itself. That is an external data source, so it needs the owner's approval under the resource policy; it has not been requested. Not a defect, an unfinished validation. |
| **P2-3** CRISPR improvements | 3.9 | **Mostly done; CFD still rejected** | Cut position (D8) and on-target exclusion (D9) are implemented and tested; the U6/Pol III poly-T consideration is named and sourced (S5). CFD scoring stays rejected: it needs the published mismatch weight matrix, and 14.3 forbids inventing it. Implementing it requires the matrix supplied as a cited data file. |
| **P2-4** cross-checks | 3.13, 10.8 | **Partially done** | REF against the reference sequence is implemented (P0-7) and the reference alphabet is now validated (S22). Primer and guide consistency against called variants is not implemented. It is feasible without external data and is the natural next functional item. |
| **P2-5** honest labels, G·U, RNAfold | 3.3, 3.4, 8.2 | **Partially done** | Unresolved residues keep their positions (Stage 5, extended by D14-F6), and the ORF label is "Heuristic score". G·U pairing and an optional RNAfold integration are not implemented; RNAfold would be an external tool and needs the dependency justification of 1.3.1. The "complexity" renaming belongs with Stage G terminology. |
| **P2-6** performance | 3.13, 11.4 | **Partially done** | The cloning planner no longer rescans (S18): 10 enzymes went from 100-110 s to 0.57 s. The mutation alignment is still quadratic in time and memory; S9 measured it at ~45 bytes per matrix cell, so 10 kb needs about 4.5 GB. A banded or linear-space alignment is a real algorithmic change and is the remaining part. |
| **P2-7** Type IIS enzymes | 3.6 "Pay particular attention to Type IIS enzymes" | **Not implemented, currently harmless** | The database holds no Type IIS enzyme, so the forward-strand search cannot produce a wrong result today. Adding one without Type IIS support would. The guard is that both arrive together. |
| **S11** live literature providers | 3.11, 1.3 | **Blocked, awaiting the owner** | Testing against live providers means sending queries to external services, which needs explicit approval. Offline fixtures cover parsing, deduplication and retries. |
| **S24** XML parsing exposure | 14.3, 1.3.1 dependency policy | **Measured, residual risk documented** | XXE is not reachable; entity expansion is capped by Expat >= 2.4.0. The cap is version-dependent, which matters because Synapse targets Python 3.10. See the audit entry. |
| **P2-2** menus inside analysis modules | 7.2 | **Partially done** | `core/api.py` now covers every analysis (S14), so a non-interactive caller never needs the menus. The menus still live in the analysis modules; moving them is a structural change of its own. |
| **P1-9** organism mode, 5'-partial, circular ORFs | 3.2 | **Deferred sub-items** | The genetic-code part is complete (RNA codon view). Organism mode, 5'-partial ORFs and ORFs spanning the origin remain unimplemented. |
| **P3-1** experiment planner | 4.3 | **Not implemented, honestly labeled** | The menu entry says "not implemented" (D15) rather than advertising a feature that does nothing. |

## Handover notes (working practice)

- Work only in this copy; never modify the original `Project Synapse`.
- 15 files carry the Windows ReadOnly attribute (`main.py`,
  `analyses/{blast,cloning,crispr,literature,mutation,orf,primer,restriction,rna,variant_analysis}.py`,
  `core/pipeline.py`, `data/{genetic_code,restriction_enzymes}.py`,
  `ui/menu.py`). Clear it only for the files a package edits and restore
  it after the commit. (`main.py` currently shows attribute 33 instead of
  524321: the OneDrive "pinned" flag is absent; only ReadOnly is managed.)
- Most source files use CRLF; `core/api.py`, the other `core/` modules
  except `pipeline.py`/`validation.py`, `fileio/genbank.py`,
  `fileio/loader.py`, the tests and the docs use LF. `analyses/orf.py` has
  exactly two intentional bare-LF lines; edit it byte-exactly.
- `pyright.exe` is blocked by a Windows application control policy
  (WinError 4551) since 2026-10-31; run it as
  `python -m pyright --pythonversion 3.10` instead. The baseline for
  comparison remains the 134 diagnostics of `ab4bd2a`; the current count is
  112 (22 of them removed with the dead primer code in D16 part 2).
- Check: `python -m unittest discover -s tests -t .`;
  `git -c core.whitespace=cr-at-eol diff --check`; pyright 1.1.414 with
  `--pythonversion 3.10` compared against the 134-diagnostic baseline of
  `ab4bd2a` (no new diagnostics allowed).
- Commit messages end with a `Co-Authored-By:` line naming the model that
  actually produced the work, not a fixed model name.
- Biopython and primer3-py are test-only; comparison tests skip without
  them. BLAST+ 2.17.0 is installed (external-tool tests skip without it);
  bcftools is not.
