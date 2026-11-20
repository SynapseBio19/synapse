# Phase 8 continuation record

**Start a new session here.** Then read:
1. `docs/PHASE8_REPORT.md`;
2. `docs/ARCHITECTURE.md`;
3. `docs/PROJECT_HISTORY_AND_DECISIONS.md` §4 (binding decisions).

## State (2026-11-18)

| Item | Value |
|---|---|
| HEAD | the commit that adds this record (after `b3ba026`) |
| Tests | 608, 0 failures on CPython 3.14.6 and 3.10.11 |
| Characterization | 182 scenarios in `tests/fixtures/characterization/`, unchanged |
| Golden master | 12/12 |
| Pyright | 11, all in `analyses/crispr.py` (allow-list) |
| Ruff | clean |
| Coverage | 83.7 % |
| ReadOnly | none (dropped by owner mandate; Git is the protection) |
| Remote | none; CI never ran on GitHub |

## Binding decisions

- **Mandate.** In the copy, Claude may:
  - restructure code;
  - fix defects that are backed by evidence;
  - delete dead code;
  - rename internals.
- **Excluded.** CRISPR and the AI Experiment Planner (which stays
  **COMING SOON**). The step a safety mechanism interrupted on 2026-11-13
  is never retried.
- **Name and domain.** The name is Synapse, the domain `synapse-gen.com`,
  the pip distribution `synapse-gen`.
- **Not chosen yet.** No licence until the owner chooses one.

## Next actions, in this order

1. Split the remaining large modules with `scripts/refactor/` (procedure
   in its README), one per commit. Order: `protein`, `variant_analysis`,
   `sequence`, `orf`, `restriction`, `rna`, `literature`.
2. Make an installable `synapse` command: package `examples/` as data,
   move `main` into a package, add `[project.scripts]`.
3. Consolidate the duplicated sequence helpers on `core.seqtools`, one
   function at a time, after recording each copy's error texts.
4. Owner steps:
   - choose the licence;
   - create the GitHub repository, so that the CI can run and its
     findings can be fixed.
