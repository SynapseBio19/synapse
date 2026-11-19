# Phase 9 continuation record

**Phase 9 finished on 2026-11-19.**
- Report: `docs/PHASE9_REPORT.md`.
- Architecture: `docs/ARCHITECTURE.md`.
- Change log: `docs/CHANGES.md`, section "Phase 9".

## State at the end

| Item | Value |
|---|---|
| Tests | 659, all passing on CPython 3.14.6 and 3.10.11 |
| Skipped on 3.14.6 | only the opt-in benchmark and the 4 opt-in refseq_rna searches |
| Coverage | 85.0 % |
| Characterization | 615 scenarios in 10 sections, bit-identical under different hash seeds |
| Golden master, Ruff | 12/12, clean |
| Pyright | 11, all in the excluded `analyses/crispr.py` |
| Remote | none; CI never ran on GitHub |

## Owner's machine

| Item | State |
|---|---|
| Smart App Control | off (since 2026-11-19) |
| BLAST databases | `C:\SynapseData\blast_db`; `SYNAPSE_BLAST_DB_ROOT` set in the user environment |
| refseq_rna | VERIFIED, release 2026-09-30; NCBI offers 2026-10-30 |
| Duplicate `data/blast_db/local_test` | in the Recycle Bin |
| Stale `refseq_rna/synapse_test.nal` | in the Recycle Bin |

## How to resume

- **Splitting a module.** Use `scripts/refactor/README.md`:
  1. characterization section;
  2. frozen names;
  3. `--check`;
  4. the split and `verify_move.py`;
  5. search the tests for patches of the public module;
  6. add the package to `pyproject.toml`, the independence test and the
     CI walk (automatic).
- **Full check.**
  - The suite on `.venv` (3.14) and `.tools/python-3.10.11`, also with
    `SYNAPSE_BLAST_DB_ROOT` set.
  - `tests.test_golden_master`, Ruff, `scripts/ci/pyright_gate.py`.
  - Optionally `SYNAPSE_REFSEQ_RNA_TESTS=1`, which takes 17 minutes.

## Open, for the owner

- **Licence.** `docs/LICENSING_AND_GOVERNANCE.md`.
- **GitHub repository.** After it exists, Claude pushes, and fixes what
  the first CI run finds.
- **Optional:** update refseq_rna to NCBI's current release with the
  database manager.
- **Excluded by owner decision:** CRISPR/CFD and the AI Experiment Planner.
