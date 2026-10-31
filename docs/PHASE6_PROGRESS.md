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
| 9 | D7 stale BLAST artifact labelled as legacy; AUDIT_STATUS cites the test | `e1a48c6` | documentation; the README's numbers match `tests/test_p0_regressions.py` | Corrected (file kept; deletion is the owner's decision) |
| 10 | D9 CRISPR on-target site excluded from the local specificity screen | `3514617` | internal tests (both strands, identical copy, near match) | Corrected; suite 203/203; pyright unchanged; S23 recorded |
| 11 | D8 CRISPR cut site per system, with citations | `ca2ec39` | internal tests (both strands, all four systems); literature (PMIDs 22745249, 25830891, 30166441, 26422227) | Corrected; suite 206/206; pyright unchanged |
| 12 | D10 dead CRISPR config fields act; D11 exports in `exports/crispr` | `8145556` | internal tests (defaults unchanged) | Corrected; suite 210/210; pyright unchanged |
| 13 | D19 pipeline `run_orf` search settings; fixed `completed` shape | `e47cb85` | internal tests (equality with `core.api.find_orfs`) | Corrected; suite 212/212; pyright unchanged |
| 14 | D15 planner entry labelled "not implemented" | `f3a09c1` | end-to-end `main.py` test | Corrected (planner itself remains P3-1, deferred); suite 213/213 |
| 15 | S15 variant-to-mutation bridge works; S16 REF_MISMATCH no longer leaks into the input | `c71963b` | internal tests; pyright 134 → 132 (the two S15 diagnostics) | Corrected; suite 215/215 |
| 16 | S17 cloning digest uses the top-strand cut | `c743beb` | internal tests; agreement with `restriction.perform_digest` for all 231 enzyme pairs, linear and circular | Corrected; suite 217/217 |
| 17 | S19 start-codon classification and scores use the genetic code | `6c20e69` | internal tests; golden master unchanged | Corrected; suite 219/219 |

| 18 | S23 PAM-proximal guide seed; S3 real RTT coverage check; S4/S5 thresholds documented and sourced | this package | internal tests; literature verified in PubMed (PMIDs 24476820, 27096365, 29160308, 31634902, 33811918) | Corrected; suite 228/228; pyright unchanged. Done by Opus 5 on a restricted assignment after a safeguard interrupted the Opus 5.5 session. |

| 19 | D16 part 1: one definition per name in `analyses/primer.py` | this package | bytecode fingerprint of all 129 module-level callables unchanged; pyright `reportRedeclaration` gone | Partially corrected; v5–v11 chain deferred; suite 229/229 |

| 20 | S21, S10 GenBank reader reports unread symbols and locations; S8 failing `to_dict()` logged; S22 variant reference must be DNA | this package | internal tests; bundled GenBank fixture unaffected | Corrected; suite 236/236 |

| 21 | S1, S2 assessed as not defects (locked by tests); S9 quantified by measurement; S11, S12 resolved in documentation | this package | execution on constructed minus-strand cases; `tracemalloc` measurement at three sizes | Assessed; suite 240/240 |

## Paused (2026-10-31)

Planned pause at the owner's request; see `docs/PHASE6_HANDOVER.md` for the
exact state, the open findings and the continuation steps. Stages F–I have
not started.

## Open after the current package

D16–D18 and S1–S6, S8–S14, S18, S21–S23 (see `docs/PHASE6_INITIAL_AUDIT.md`)
until a package above closes an item. D12, D1–D11, D13–D15, D19, S7, S15–S17, S19
and S20
are corrected.

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
  131 (three corrected, none added).
- Check: `python -m unittest discover -s tests -t .`;
  `git -c core.whitespace=cr-at-eol diff --check`; pyright 1.1.414 with
  `--pythonversion 3.10` compared against the 134-diagnostic baseline of
  `ab4bd2a` (no new diagnostics allowed; 132 since the S15 fix).
- Commit messages end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Biopython and primer3-py are test-only; comparison tests skip without
  them. BLAST+ 2.17.0 is installed (external-tool tests skip without it);
  bcftools is not.
