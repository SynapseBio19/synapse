# Phase 7 continuation record

Phase 7, engineering stabilization, under the owner's Phase 7 Master
Execution Prompt (2026-11-17). This file records where the work stands
so that a new session continues from it instead of starting over.

## Binding decisions

- The AI Experiment Planner stays **COMING SOON**.
- CRISPR- and CFD-related work is excluded: no development, research,
  licensing questions or retry of the step a safety mechanism interrupted
  on 2026-11-13. `analyses/crispr.py` is not edited; its pyright findings
  are left and reported.
- No new scientific features. Public APIs, workflows and results are
  preserved; any behavior change is justified, tested and documented.
- Licence: no LICENSE file until the owner chooses one.
- Safeguards: work only in this copy; ReadOnly on the 15 protected
  files; no history rewrite; Windows application control is not bypassed.

## Sources inspected for the history

- The repository records: `docs/PROJECT_HISTORY_AND_DECISIONS.md`, the
  Phase 6 final report, final assessment, continuation record,
  FEATURE_STATUS, AUDIT_STATUS, CHANGES, README and METHODS. The Phase 6
  Master Prompt is in the owner's Downloads folder.
- Claude Code sessions:
  - "Project Synapse Phase 6" and "Project Synapse Handover Setup": read
    in full for owner messages on 2026-11-13; no new activity since;
  - the current session, 2026-11-13 to 2026-11-17.
- Not accessible: the claude.ai Project conversations and the earlier
  "Claude Cowork" cycle.

## Baseline (2026-11-17, HEAD `cf930cf`)

| Item | Value |
|---|---|
| Branch / tree | `main`, clean, 155 commits, no Git remote |
| Tests | 516, 0 failures, 5 skipped (opt-in benchmark; four primer3-py comparisons because Windows application control blocks its compiled module) |
| Golden master | 12/12 |
| Pyright 1.1.414 (`python -m pyright --pythonversion 3.10`) | 58 errors |
| Pyright errors, by location | tests 27, `analyses/crispr.py` 11 (excluded), `literature.py` 8, `rna.py` 4, `orf.py` 3, `protein.py` 2, `core/report.py` 2, `blast.py` 1 |
| Pyright errors, by rule | `reportArgumentType` 19, `reportOptionalMemberAccess` 18, `reportReturnType` 5, `reportPossiblyUnbound` 3, `reportAttributeAccessIssue` 3, `reportCallIssue` 3, `reportOptionalOperand` 3, `reportOperatorIssue` 2, `reportMissingImports` 1, `reportIncompatibleMethodOverride` 1 |
| Ruff 0.14.4 | pyflakes (F) clean; with the default rule set 20 findings (E701 15, E741 3, E402 1, E731 1); no Ruff configuration in `pyproject.toml` |
| Coverage | 75 % combined (coverage.py 7.10.7) |
| Python | 3.14.6 (`.venv`), 3.10.11 (`.tools/`, embeddable); no other interpreters |
| Packages (`.venv`) | biopython 1.88, coverage 7.10.7, numpy 2.5.4, primer3-py 2.4.0, pyright 1.1.414, ruff 0.14.4, ViennaRNA 2.8.0 |
| ReadOnly | 15/15 |

## Progress

| Package | State | Commit |
|---|---|---|
| A Documentation accuracy | done | (this commit) |
| B Pyright | open | |
| C Refactoring | open | |
| D Tests and coverage | open | |
| E CI | open | |
| F Dependencies and installation | open | |
| G Error handling | open | |
| H Reproducibility | open | |
| I Licence and governance | open | |
| J Final validation and report | open | |

## Next action

Package B: triage the pyright findings outside `analyses/crispr.py`,
starting with `analyses/literature.py` and `core/report.py`.
