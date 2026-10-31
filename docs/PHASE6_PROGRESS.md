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
| 2 | Central alphabets and classifier (F12/D13 groundwork) | this package | internal tests; external reference comparison (`Bio.Data.IUPACData`); literature (PMIDs 2582368, 6692818) | Added, not yet wired in; suite 167/167; pyright unchanged |

## Open after the current package

See `docs/PHASE6_INITIAL_AUDIT.md` (D6–D11, D13–D19, S1–S22) until a
package above closes an item.
