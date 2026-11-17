# Phase 7 continuation record

Phase 7, engineering stabilization, under the owner's Phase 7 Master
Execution Prompt (2026-11-17). **Finished on 2026-11-17.** Final report:
`docs/PHASE7_FINAL_REPORT.md`. Verdict: PARTIALLY COMPLETED.

- **Done:** documentation, type safety, CI configuration, dependencies,
  installation, messages and traceability.
- **Partial:** refactoring (only `main()`) and coverage (81 % of the 85 %
  target).
- **Blocked or deferred:** CI on a server (no remote), licence, the
  legacy files and unreachable printers (owner decisions), and
  CRISPR/CFD and the AI Experiment Planner (owner decisions).

## Binding decisions

- The AI Experiment Planner stays **COMING SOON**.
- CRISPR- and CFD-related work is excluded. `analyses/crispr.py`,
  `ui/crispr_menu.py` and the CRISPR tests are not edited; the 11 pyright
  findings in `crispr.py` are on the allow-list of
  `scripts/ci/pyright_gate.py`. The step a safety mechanism interrupted on
  2026-11-13 is not to be retried.
- No new scientific features. Public APIs, workflows and results are
  preserved.
- No LICENSE file until the owner chooses
  (`docs/LICENSING_AND_GOVERNANCE.md`).
- Safeguards:
  - work only in this copy;
  - ReadOnly on the 15 protected files;
  - no history rewrite;
  - Windows application control is not bypassed (it blocks primer3-py's
    DLL on this machine).

## State at the end (2026-11-17)

| Item | Value |
|---|---|
| HEAD | the commit that adds the final report |
| Working tree | clean |
| Commits | 167 |
| Tests | 565, 0 failures on CPython 3.14.6 (5 skipped) and 3.10.11 (23 skipped) |
| Golden master | 12/12 |
| Pyright | 11 errors, all in `analyses/crispr.py` |
| Ruff | clean (configured in `pyproject.toml`) |
| Coverage | 81 % |
| ReadOnly | 15/15 |

## Commits of the phase

| Commit | Package | Content |
|---|---|---|
| `96d0f07` | A | FEATURE_STATUS checked row by row; 8 overclaimed functions tested; documentation-consistency test |
| `2f8f6e6` | B | Pyright 58 -> 11 |
| `cdbb7ac` | D/G | BLAST dialog tests; `.njs` left behind by Remove Database fixed |
| `bedcd42` | D | `main()` driven in-process; single-primer and overview tests |
| `ca5d53a` | D/G | Literature pipeline tests; BibTeX escaping fixed |
| `c97695d` | C | Ruff default rules clean, Ruff configured |
| `895fc5e` | C | `main()` split into session steps, output byte-identical |
| `37986d7` | E/F | GitHub Actions workflow, gate scripts, `dev` extra, clean install |
| `3f48f1f` | G | Message prefixes and rules |
| `22d0abb` | H | Genetic code and input mode in RNA and protein report entries |
| `236008f` | I | Licence recommendation, CONTRIBUTING |
| (final) | J | Final verification and reports |

## Next actions (only on the owner's request)

1. Owner decisions:
   - the licence and copyright holder;
   - keep or delete the two legacy BLAST result files and the two
     unreachable report printers.
2. Push the repository to a hosted remote so that the CI runs. Then fix
   whatever it finds on Linux, Python 3.11-3.13 and the reference job
   (primer3-py runs there).
3. Split the large modules in the order given in the final report,
   section E.
4. Raise coverage in the restriction, cloning and protein dialogs.
