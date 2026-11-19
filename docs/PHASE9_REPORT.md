# Phase 9 report: completion under the owner's mandate

Dates: 2026-11-18/19.
- **Mandate.** Fix every defect and complete the project, so that only the
  organizational steps (GitHub, licence) remain for the owner.
- **Excluded.** CRISPR and the AI Experiment Planner.
- **Precondition.** Windows Smart App Control, which had blocked primer3-py
  and ViennaRNA, was switched off by the owner on 2026-11-19.

## 1. Summary

Phase 9 found and fixed seven defects:
- three in the science (BLAST E = 0, the hash-seed-dependent sequence
  complexity, missing genetic codes);
- four in the product (runtime views left behind, user files in
  site-packages, a self-test that failed with any configured database
  folder, a CI quoting error).

It added:
- an installable `synapse` command;
- a size, release and disk check before database downloads;
- searches against the owner's real refseq_rna database.

All large analysis modules except `literature` now have a structure by
responsibility. Every split was frozen first in a characterization
snapshot. Test coverage reached the 85 % target.

| Metric | Start of Phase 8 (`9fad47e`) | End of Phase 9 |
|---|---:|---:|
| Tests | 601 | 659, all passing |
| Coverage | 82.3 % | 85.0 % |
| Characterization scenarios | 78 | 615 in 10 sections |
| Genetic codes | 4 | 27 (all of NCBI) |
| Largest analysis module (not excluded) | 4,670 lines | 2,095 (`literature`, kept whole on purpose) |
| Pyright | 11 (crispr.py only) | 11 (crispr.py only) |
| Ruff, golden master | clean, 12/12 | clean, 12/12 |

## 2. Defects found and fixed

| # | Defect | How it was found | Evidence of the fix |
|---|---|---|---|
| 1 | BLAST: an E-value of 0 was read as infinity. The best hits were rated NOT_SIGNIFICANT, ranked last and exported as invalid JSON (local and remote searches) | first search against the real refseq_rna | `ZeroEvalueTests`; the real search now gives E = 0 and VERY_STRONG |
| 2 | Local BLAST left its runtime views (hard links) behind; the database list showed a phantom 81-GB copy | the owner's database folder | `RuntimeViewCleanupTests`; no view left after 4 real searches |
| 3 | Sequence complexity depended on Python's hash seed; the reported complexity maximum moved between runs | recording the characterization snapshot | `tests/test_determinism.py`; all 10 sections bit-identical under 3 seeds |
| 4 | An installed Synapse would write logs, exports and BLAST databases into site-packages | building the `synapse` command | `StorageLocationTests`; clean install writes to the user folder |
| 5 | The BLAST self-test failed whenever `SYNAPSE_BLAST_DB_ROOT` pointed outside `data/blast_db`, which was the owner's setting | running the suite with the owner's environment | the self-test checks the exact root per configuration; the suite runs in a fixed environment |
| 6 | The CI import check had a shell-quoting error (introduced by me in Phase 8) | running the workflow lines exactly as written | the step now passes on a clean install |
| 7 | The session ended with a Python traceback on Ctrl+C or end of input | the command test | `cli()`; `tests/test_app_entry.py` |

## 3. Completed features

- **Genetic codes.** All 27 NCBI codes, checked against NCBI's `gc.prt`
  4.6 (stored with its checksum) and against Biopython. The
  context-dependent codons of codes 27, 28 and 31 are handled and
  documented.
- **BLAST databases.**
  - Before a download, Synapse shows NCBI's release date, size, the space
    needed and the free space. It refuses the download when the space is
    too small, and says when NCBI's release is newer than the installed
    one.
  - The owner's refseq_rna was moved out of OneDrive to
    `C:\SynapseData\blast_db`.
- **Installation.**
  - `pip install .`, then `synapse`; also `python -m synapse_gen`, and
    `python main.py` from the source folder.
  - The example file is installed with the package.
  - User files go to a per-user folder (`SYNAPSE_HOME`).
  - The CI imports all 138 installed modules.

## 4. Verified against references

- **refseq_rna, real searches.**
  - Exact, three substitutions (99.25 %), reverse complement, and a 22-nt
    query found by the rescue: `tests/test_blast_refseq_rna.py`, opt-in,
    passed in 17.5 minutes.
  - A live NCBI search of the same query was still queued at NCBI after
    25 minutes (NCBI's own status page: "Searching"). Synapse read the
    status correctly and offered the RID for resuming.
  - The remote result therefore could not be compared. Remote searches use
    the same XML parser as the verified local ones.
- **ORF finder.** Agrees with an independent naive finder on 675 random
  sequences in all 27 codes.
- **Biopython comparisons.** No defect found:
  - the seven reverse-complement implementations;
  - the optimality of the mutation alignment;
  - GC content and RNA translation in all codes.
- **primer3-py and ViennaRNA comparisons** run again (Smart App Control
  off).

## 5. Structure

The public modules stay the interface. The code is in internal packages,
each split frozen first and moved verbatim (`docs/ARCHITECTURE.md`,
`scripts/refactor/`):
- `protein`;
- `variant_analysis`;
- `sequence`;
- `orf`;
- `restriction`;
- `rna`.

**Kept whole, on purpose.** `analyses/literature.py`: its tests replace
module functions 20 times and the providers share a throttle, so a split
would let those tests pass without reaching the code.

## 6. What remains, and for whom

| Item | Who | Note |
|---|---|---|
| Licence | owner | options in `docs/LICENSING_AND_GOVERNANCE.md`; not chosen on the owner's behalf |
| GitHub repository, first CI run | owner creates the repository; Claude can then push and fix what the CI finds on Linux and Python 3.11-3.13 | the workflow is ready and was run step by step locally |
| Update of refseq_rna to NCBI's release of 2026-10-30 | owner's choice | the database manager now offers it with size and space |
| Scientific benchmarks on curated data, laboratory validation | owner / laboratory | outside software verification |
| CRISPR (CFD, genome-wide off-targets), AI Experiment Planner | owner decision | excluded |
| Graphical or web interface | owner decision | not part of the mandate |

## 7. Assessment

The software is now in a state where:
- every analysis module is structured by responsibility and protected
  against silent changes, through exact snapshots, namespace and
  independence tests and a determinism test;
- seven real defects found by testing against real data and real
  environments are fixed;
- the project installs and runs as a command, keeps user data in the
  right place and checks resources before large downloads.

The remaining steps are organizational (licence, repository), scientific
(benchmarks on curated data, laboratory work) or excluded by the owner.
