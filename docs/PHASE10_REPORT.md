# Phase 10 report: genome scale, curated benchmarks, release preparation

**Date:** 2026-11-20. **Start:** `a2898aa` (end of Phase 9).
**Owner's request:** do everything that can still be done in the project,
as professionally as possible, and report when Synapse is finished on the
AI side; the licence and GitHub follow from the owner. The mandate of
2026-11-18 holds; CRISPR/CFD and the AI Experiment Planner stay excluded.

## A. Summary

Synapse now runs every analysis on a complete bacterial genome, has been
benchmarked against curated reference data (RefSeq gene models, ClinVar
variants), has a user guide, and is prepared for its first public release.
Seven defects were found and fixed, each with a test that fails on the
old code. Nothing in the excluded areas was touched.

## B. Work packages

| Package | Commit | Content |
|---|---|---|
| A | `fb5f978` | Genome-length input: ORF relationships, restriction annotation and site search, sequence and RNA window scans, hairpin screen, profile thinning above 200,000 windows |
| B | `28d4db9` | HGVS benchmark on 1,500 ClinVar variants; allele 3' shift fix |
| C | `870c91c` | Per-ORF features computed once |
| D | `b3f9f5c` | Dialog tests for protein, restriction, cloning; two dialog fixes |
| E | `51fbba4` | Literature author de-duplication across PubMed and Europe PMC |
| F | `f70e665`, `66554aa` | Literature snapshot, then the split into `analyses/_literature/`; splitting-tool fix |
| G | `cf30674` | User guide, pUC19 example, release notes, `CITATION.cff`, GitHub templates |
| H | this commit | Benchmarks document, method and status documentation, this report |

## C. Defects found and fixed

| # | Area | Defect | Found by | Evidence |
|---|---|---|---|---|
| 1 | ORF | All-pairs overlap analysis: 8 minutes on *E. coli* | genome benchmark | `test_orf_relationship_sweep` |
| 2 | Restriction | Fragment annotation scanned every cut per fragment: no result in 15 minutes | genome benchmark | `test_restriction_fragment_index` |
| 3 | Sequence overview | 3.5 GB for the profiles of a genome; minutes per scan | genome benchmark | `test_window_scans` |
| 4 | RNA | Hairpin screen and profiles over 10 minutes on a genome | genome benchmark | `test_rna_window_scans` |
| 5 | HGVS | Multi-change alleles could overlap after the 3' shift (PTEN c.528_534delinsA) | ClinVar benchmark | `test_hgvs_allele_shift` |
| 6 | Literature | Authors listed twice when PubMed and Europe PMC records merge | literature snapshot | `test_literature_authors` |
| 7 | UI | Three protein errors without the ERROR prefix; Back inside a cloning analysis left the cloning menu | dialog tests | `test_protein_restriction_cloning_dialogs` |

Also fixed: the splitting tool generated only one of several
`import urllib.x` lines (found by Pyright).

Defects 1-4 are performance defects: the results were correct but not
obtainable for genome-length input. The faster code returns the same
results, bit for bit (section D).

## D. Evidence that results did not change

- Differential tests against the former implementations (five new test
  modules).
- End to end: seven analyses on *E. coli* segments give identical SHA-256
  digests of their complete output with the code of `a2898aa` and with
  the current code (docs/BENCHMARKS.md).
- Characterization: 636 scenarios in eleven sections, unchanged; golden
  master 12/12 unchanged.
- The one intended change of output is documented: above 200,000 profile
  windows the stored profiles are thinned (`profile_resolution`).

## E. Benchmarks

Details: `docs/BENCHMARKS.md`.

- Genome: overview 154 s → 41 s (3.5 GB → 0.64 GB); ORF search
  499 s → 12 s; RNA > 600 s → 137 s; digest with 22 enzymes > 900 s → 6 s.
- ORFs vs RefSeq: 4,250 of 4,294 CDS at their stop codon; every CDS of at
  least 90 nt.
- HGVS vs ClinVar: c. 1,489/1,500 identical, 11 equivalent; p. 1,186/1,186.
- Not done: Tm against measured values (no licensed dataset obtained).

## F. Tests and quality

| Check | Result |
|---|---|
| Suite, CPython 3.14.6 | 722 tests, 0 failures; 5 opt-in skips (refseq_rna searches, 10-kb benchmark) |
| Suite, CPython 3.10.11 | 722 tests, 0 failures; 29 skips (the opt-in tests and the optional packages not installed for 3.10) |
| Golden master | 12/12 unchanged |
| Characterization | 636 scenarios in 11 sections, unchanged |
| Ruff | clean |
| Pyright | 11 errors, all in the excluded `analyses/crispr.py` |
| Coverage (line + branch, combined) | 85.0 % → 86.1 %; cloning dialog 74 → 82 %, restriction 75 → 85 %, protein 83 → 87 % |
| CI import step, exactly as in `ci.yml`, on a clean non-editable install | 150 modules import from site-packages; self-tests pass; `synapse` starts |
| Full ClinVar benchmark (`SYNAPSE_CLINVAR_FULL=1`) | 5 tests pass on all 1,500 variants |
| Live NCBI BLAST against `refseq_rna` | not completed: NCBI did not finish within 45 minutes (RID FZ5HPRV0014); open, see docs/BENCHMARKS.md |

New test modules: `test_window_scans`, `test_rna_window_scans`,
`test_orf_relationship_sweep`, `test_orf_feature_speed`,
`test_restriction_fragment_index`, `test_hgvs_allele_shift`,
`test_hgvs_clinvar_benchmark`, `test_protein_restriction_cloning_dialogs`,
`test_literature_authors`, `test_user_guide_examples` (659 → 722 tests).

## G. Documentation

New: `docs/USER_GUIDE.md`, `docs/BENCHMARKS.md`, `docs/RELEASE_NOTES.md`,
`docs/PHASE10_REPORT.md`, `CITATION.cff`, `.github/ISSUE_TEMPLATE/*`,
`.github/PULL_REQUEST_TEMPLATE.md`, `.github/dependabot.yml`,
`tests/fixtures/clinvar/SOURCE.md`.
Updated: `README.md`, `docs/METHODS.md`, `docs/FEATURE_STATUS.md`,
`docs/CHANGES.md`, `docs/ARCHITECTURE.md`, `docs/TOOL_AND_DATA_LEDGER.md`,
`docs/TEST_COVERAGE.md`, `docs/PROJECT_HISTORY_AND_DECISIONS.md`,
`scripts/refactor/README.md`.

## H. What is left, and whose it is

**The owner's** (`docs/RELEASE_NOTES.md`, "Before the first release"):
1. Choose the licence; add `LICENSE`, the `license` fields in
   `pyproject.toml` and `CITATION.cff`.
2. Create the GitHub repository and push. The CI has never run on a
   server; its first run on Linux may find platform differences, which
   Claude can fix once the repository exists.
3. Version 0.3.0, tag, release; optionally Zenodo and PyPI.
4. A contact address, if a security policy or code of conduct is wanted.

**Excluded by owner decision:** CRISPR/CFD (no CFD score, no genome-wide
off-targets), the AI Experiment Planner (Coming Soon), a graphical
interface.

**Open checks that depend on an external service:**
- The live NCBI BLAST search against `refseq_rna` (queued over 45 minutes
  at NCBI on 2026-11-20, as on 2026-11-19).

**Open scientific points, each needing data or a decision:**
- Primer Tm against measured melting temperatures (a licensed dataset).
- The instability weight KQ (24.64 here, 24.68 at ExPASy; the 1990 paper).
- Affine gap costs for the mutation alignment, which would write the 11
  split ClinVar delins as single changes (a model change).
- Experimental validation of any result (a laboratory).

## I. Assessment

On the AI side, the work that can be done without the owner's decisions,
external data or a laboratory is done. Synapse is a tested research and
teaching workbench: its calculations agree with established tools and
curated data wherever they were compared, it handles genome-length input,
it says what each result means and does not mean, and it is documented
for users and contributors. It is not experimentally or clinically
validated, and it does not claim to be.
