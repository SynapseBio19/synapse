# Phase 7.1 continuation record

Phase 7.1, architectural modularization, under the owner's Phase 7.1
prompt (2026-11-18). **Finished on 2026-11-18.** Final report:
`docs/PHASE71_FINAL_REPORT.md`; architecture: `docs/PHASE71_ARCHITECTURE.md`.

## State at the end

| Item | Value |
|---|---|
| HEAD | the commit that adds this record (after `f629add`) |
| Working tree | clean (ignored: `build/`, `synapse.egg-info/` from the install check) |
| Tests | 601, 0 failures on CPython 3.14.6 (1 skipped) and 3.10.11 (23 skipped) |
| Golden master | 12/12 |
| Characterization snapshot | 78 scenarios, unchanged since `af7cd0a` |
| Pyright | 11 errors, all in `analyses/crispr.py` (allow-list) |
| Ruff | clean |
| Coverage | 82.3 % (BLAST 84.9 %, primer 91.7 %) |
| ReadOnly | at the end of Phase 7.1: the 15 original files plus the 27 files of the internal packages. **Dropped on 2026-11-18 (Phase 8, owner mandate); Git is the protection.** |

## Completed stages

1. Baseline verified.
2. Architecture audit.
3. BLAST split.
4. Primer: explicit configuration, then the split.
5. Quality: component tests, coverage, the CI import check.
6. Documentation.
7. Final validation.

## How to work with the new layout

- **Imports.** Import from `analyses.blast` and `analyses.primer` only.
- **Changing code.** Change it in the internal module named in
  `docs/PHASE71_ARCHITECTURE.md` §3/§4. (The ReadOnly flags were dropped
  in Phase 8.)
- **Patching in tests.** Patch the module that calls the function, for
  example `analyses._blast.ncbi.urlopen`.
- **A new public name** must be added to the internal module and to the
  re-export list of the public module. `test_module_namespaces` guards the
  existing names only.
- **A new internal package** must be added to `packages` in
  `pyproject.toml`. The CI import check fails without it.
- **Results.** `tests/test_characterization.py` must stay green. Rewrite
  its snapshot (`SYNAPSE_UPDATE_CHARACTERIZATION=1`) only with a reason in
  `docs/CHANGES.md`.

## Binding decisions (unchanged)

- The AI Experiment Planner stays **COMING SOON**.
- CRISPR/CFD work is excluded. `analyses/crispr.py`, `ui/crispr_menu.py`
  and the CRISPR tests are not edited.
- No licence until the owner chooses (`docs/LICENSING_AND_GOVERNANCE.md`).
- No history rewrite, no remote, no push without the owner.

## Next actions (only on the owner's request)

1. Owner decisions:
   - licence;
   - the two legacy BLAST files and the two unreachable printers;
   - whether ReadOnly should keep covering the internal packages.
2. Push to a hosted remote so that CI runs.
3. Apply the same split to `analyses/mutation.py`, then
   `analyses/cloning.py`:
   - namespace test and characterization snapshot first;
   - then the move;
   - then the independence test.
4. Consolidate the duplicated sequence helpers on `core.seqtools`, one
   function at a time (audit F6).
