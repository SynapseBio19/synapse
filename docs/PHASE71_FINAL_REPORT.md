# Phase 7.1 final report: architectural modularization

Date: 2026-11-18. Governed by the owner's Phase 7.1 prompt. Start:
`af7cd0a`. The commits are listed in section L.

## A. Executive summary

- **BLAST and primer modularized.** `analyses/blast.py` (4,667 lines)
  and `analyses/primer.py` (5,182 lines) now keep their code in the
  internal packages `analyses/_blast/` (14 modules) and
  `analyses/_primer/` (11 modules). Each module has one responsibility,
  and the import graph has no cycle. The two public modules remain the
  only interface; they re-export every name they had.
- **Code moved verbatim.** A script cut the modules at statement
  boundaries. An independent check found all 274 statements unchanged in
  their new place. The only edit inside moved code concerns four path
  lookups that now name the `analyses` directory explicitly; the paths
  themselves are unchanged.
- **Results unchanged.** A characterization snapshot of 78 scenarios was
  recorded from `af7cd0a` before any move and is compared exactly. The
  golden master, the external reference comparisons and the
  byte-identical rendered menus are unchanged too.
- **A hidden global removed.** The primer screen read its configuration
  from a module global that each design overwrote. It is now an argument.
- **Three latent defects found by the work.**
  1. `pyproject.toml` lists packages explicitly. Without an update,
     `pip install .` would have shipped a `blast.py` that cannot import.
  2. CI ran its import check in the checkout, so it could not have
     noticed (1). It now checks the installed package.
  3. The rule tests (layering, UI rules, spelling, vocabulary, feature
     inventory) scanned `analyses/*.py` only and would have silently
     exempted moved code. They now include `analyses/**`.
- **Tests and coverage.** 565 -> 601 tests. Coverage 81 % -> 82.3 %;
  BLAST 77 % -> 84.9 %, primer 86 % -> 91.7 %.

**What materially improved.**
- A change to the BLAST XML parser, the NCBI client, the local database
  manager or the primer thermodynamics is now made, read and reviewed in
  a module of 100-1,600 lines that holds only that responsibility, not in
  a 5,000-line file.
- Coverage is visible per responsibility. That is how the
  least-tested parts were found and tested.

## B. Verified baseline (before any change, 2026-11-18)

| Item | Value |
|---|---|
| Git | `main`, HEAD `af7cd0a`, 167 commits, clean, no remote |
| Python | 3.14.6 (`.venv`), 3.10.11 (`.tools/python-3.10.11`) |
| Tests | 565, 0 failures; 3.14.6: 1 skipped; 3.10.11: 23 skipped |
| primer3-py | loads again; its four comparisons **ran and passed** (they skipped on 2026-11-17 because Windows application control blocked the module) |
| Coverage | 81 % |
| Pyright | 11 errors, all in `analyses/crispr.py` (allow-list) |
| Ruff | clean |
| Golden master | 12/12 |
| ReadOnly | 15/15 |

The Phase 7 report's figures were confirmed. Previous Claude conversations
were not re-read in this session; the project history is taken from
`docs/PROJECT_HISTORY_AND_DECISIONS.md`, the Phase 7 records and this
conversation's summary.

## C. Architecture audit (`docs/PHASE71_ARCHITECTURE_AUDIT.md`)

Main findings, by priority:

| ID | Priority | Finding | Outcome |
|---|---|---|---|
| F1 | P1 | `blast.py`: ten responsibilities in one module, a sixth of it self-test | done (D) |
| F2 | P1 | `primer.py`: historical layers, calls 3,000 lines forward | done (E) |
| F3 | P1 | hidden global configuration in the primer screen | done (E) |
| F4 | P2 | self-tests inside production modules | done for BLAST and primer; others open |
| F5 | P2 | `mutation.py`, `cloning.py` not examined in depth | open |
| F6 | P2 | duplicated sequence helpers (7 x `reverse_complement`) | open |
| F7 | P2 | dialog functions of 400-570 lines | open |
| F8-F11 | P3 | unreachable printers, import-time configuration, tests coupled to module layout, private helpers used by the UI | recorded; F10 fixed for the tests touched |

Risks identified before moving, each handled:
- the `BlastRun` and `rank_hits` cycle;
- `Path(__file__)` lookups;
- tests patching module attributes;
- primer forward references;
- `lru_cache` identity;
- the hidden global.

## D. BLAST modularization

- **Files.** `analyses/blast.py` became the public module (217 lines). New:
  `analyses/_blast/` with `__init__` plus `paths`, `config`, `models`,
  `query`, `parsing`, `ncbi`, `executables`, `runtime`, `databases`,
  `local_engine`, `evidence`, `reporting`, `workflow`, `selftest`.
  Responsibilities and dependencies: `docs/PHASE71_ARCHITECTURE.md` §3.
- **Extracted.**
  - Configuration and catalogs; the data model; query validation; XML
    parsing.
  - The NCBI client; executable discovery; runtime alias views; the
    database manager; local execution.
  - Evidence and ranking; serialization, exports and reports; the
    workflow; the self-test.
- **Compatibility.**
  - All 125 names are reachable (`test_module_namespaces`).
  - The public names are the implementation objects.
  - Each internal module imports on its own.
  - The snapshot (run, report text, CSV/JSON/text exports, evidence,
    comparison, empty database manager, database root, executable search
    directories, user-path base, self-test) is unchanged.
  - Rendered menus are byte-identical. A clean non-editable install works.
- **Changed for compatibility.** Two tests patched names that
  `analyses.blast` only imported (`time.sleep`, `urlopen`); they now patch
  `analyses._blast.ncbi`. Such names (standard library, `core.alphabets`)
  are no longer attributes of `analyses.blast`.
- **Tests.**
  - The existing BLAST tests (11 modules), the snapshot and the namespace
    tests.
  - New `test_blast_components` (19): validation messages, the BLAST+
    search order, local-engine failure paths and the short-query rescue.

## E. Primer modularization

- **Files.** `analyses/primer.py` became the public module (244 lines,
  both `__main__` blocks). New: `analyses/_primer/` with `__init__` plus
  `sequence`, `models`, `thermodynamics`, `structure`, `specificity`,
  `design`, `professional`, `screening`, `ranking`, `workflow`,
  `selftest`. Details: `docs/PHASE71_ARCHITECTURE.md` §4.
- **Extracted.**
  - Sequence helpers; the legacy model and engine; thermodynamics (all
    Tm methods, parameter sets, salt terms, presets).
  - Structure screens; the specificity screen.
  - The professional model and scores; the candidate screen; pair
    ranking; the professional workflow and restriction-site design.
  - The self-tests.
- **Done first, in its own commit.** The configuration of the screen is
  passed explicitly, and `_V5_ACTIVE_CONFIG` is removed.
- **Compatibility.**
  - All 139 remaining names are reachable.
  - `lru_cache` wrappers remain single objects.
  - The one statement that defines no name (the 2004 AA/TT stack update)
    sits in `thermodynamics` right after its dict. A first run of the
    script had left it in the public module; the new independence test
    covers this.
  - The snapshot is unchanged: legacy design, professional design under
    three Tm models and high salt, restriction-site design, single-primer
    analyses, Tm comparisons, structure checks, condition summaries, four
    self-tests.
- **Primer3.** The primer3-py comparisons (Tm across salt conditions,
  structure ordering) **ran and passed** on 2026-11-18, before and after
  the split.
- **Tests.**
  - The existing primer suites, including the screen-equivalence tests.
  - New `test_primer_screen_configuration` (3) and
    `test_primer_components` (8): condition validation, the specificity
    screen, workflow warnings and failure reasons.

## F. Testing

| Metric | Before (`af7cd0a`) | After |
|---|---:|---:|
| Tests | 565 | 601 |
| Failures | 0 | 0 |
| Skipped, 3.14.6 / 3.10.11 | 1 / 23 | 1 / 23 |
| Coverage, total | 81 % | 82.3 % |
| BLAST / primer coverage | 77 % / 86 % | 84.9 % / 91.7 % |
| Golden master | 12/12 | 12/12 |

**New tests (36).**
- `test_module_namespaces` (4)
- `test_characterization` (2 tests, 78 scenarios)
- `test_primer_screen_configuration` (3)
- `test_blast_components` (19)
- `test_primer_components` (8)

**Regressions.** Found and fixed before any commit: the statement left in
the public module (E), the packaging list (A.1) and the narrow rule scans
(A.3). No test failed in a committed state.

**Not reached.** 85 %. The rest of the gap is in `ui` (74 %): the CRISPR
dialog is excluded, `print_protein_report` is unreachable pending the
owner's decision, and BLAST download and remote paths remain.

**Observation recorded, not changed.** On a poly-A template the
professional design returns one A21/T21 pair, flagged high and critical
risk. This is the documented "measure everything, flag in the ranking"
design.

## G. Code quality

- **Pyright.** 11 errors before and after, all in the excluded
  `analyses/crispr.py`. No new error; the new tests type-check.
- **Ruff.** Clean. `analyses/blast.py` and `analyses/primer.py` ignore
  F401 per file, because their imports are the interface.
- **Type safety.** Unchanged in the moved code (moved verbatim). The
  explicit screen configuration removes the one piece of hidden state in
  the primer engine.
- **Python 3.10.** The full suite passes on 3.10.11.
- **Line endings.** The public modules keep CRLF; the new files use LF,
  like all files added since Phase 6.

## H. CI

- **Change.** The import check of the `tests` job runs in `$RUNNER_TEMP`
  with `python -I`. It asserts that the modules come from `site-packages`
  and runs the BLAST self-test there.
- **Checked locally.**
  - The commands were run against a clean non-editable install: they
    pass.
  - With `analyses/_primer` removed from that installation, the step fails
    as intended.
  - Every other job's commands ran locally (suite, golden master, Ruff,
    pyright gate, coverage).
  - The YAML was reviewed by hand; no YAML parser is installed, and none
    was added.
- **Remote.** **Unverified.** The repository has no remote; nothing was
    pushed, no repository was created. A local run is not a GitHub
    Actions run.

## I. Documentation

**Created.**
- `docs/PHASE71_ARCHITECTURE_AUDIT.md`
- `docs/PHASE71_ARCHITECTURE.md`
- `docs/PHASE71_FINAL_REPORT.md`
- `docs/PHASE71_CONTINUATION.md`

**Updated.**
- `docs/CHANGES.md` (Phase 7.1 section)
- `docs/FEATURE_STATUS.md` (test count, primer3-py row)
- `docs/AUDIT_STATUS.md`
- `docs/TEST_COVERAGE.md`
- `docs/PROJECT_HISTORY_AND_DECISIONS.md` (§10)
- `README.md` (structure)
- comments in `.github/workflows/ci.yml` and `pyproject.toml`

## J. Remaining technical debt

- **Unfinished.**
  - `mutation.py`, `cloning.py`, `protein.py`, `variant_analysis.py` and
    `sequence.py` are not split (F5); the same method applies.
  - Duplicated helpers (F6) and long dialogs (F7) are open.
- **Architectural limitations.**
  - `LocalBlastDatabaseManager` (985 lines) and `_primer/thermodynamics`
    (1,581) are cohesive but long.
  - BLAST configuration is read at import time (F9).
  - The `_v5`-`_v11` names are kept, because they are part of the tested
    private surface.
- **Risks.**
  - The internal packages are new to the owner's workflow: ReadOnly now
    also covers their 27 files.
  - Tests must patch the implementing module (documented in
    `PHASE71_ARCHITECTURE.md` §2).
  - CI has never run remotely.
- **Housekeeping.** The local install checks recreated the ignored `build/`
  and `synapse.egg-info/` directories. They were deleted after the first
  check and left in place after the second. Both are build artifacts,
  covered by `.gitignore`.

## K. Explicitly deferred areas

- **CRISPR/CFD.**
  - Not modified, refactored, tested, investigated or extended.
  - `analyses/crispr.py` keeps its 11 allow-listed pyright findings.
  - CRISPR depends on neither BLAST nor primer.
  - The rule-scan helpers in `tests/test_phase6_tier2.py` were widened;
    the CRISPR tests in that file are unchanged.
- **AI Experiment Planner.** Still **COMING SOON**, unchanged.
- **Other.**
  - No licence chosen.
  - No new scientific feature; no algorithm changed.
  - No UI change: rendered menus are byte-identical.
  - No new dependency.
  - No remote.

## L. Final assessment

**Engineering maturity after Phase 7.1.** The two largest modules have a
structure that matches their responsibilities. It is protected by tests
that check names, exact results and module independence, and the moves
were shown to be verbatim. Packaging and CI have one blind spot fewer.

This is a substantive improvement in maintainability. It is not a
claim of production readiness: the remaining large modules, the open
duplication, 74 % UI coverage and a CI that has never run on a server
are real.

| Status | Items |
|---|---|
| Implemented and tested | BLAST and primer modularization; explicit screen configuration; namespace, independence and characterization tests; component tests; packaging fix; widened rule scans; primer3-py comparisons run |
| Implemented, not fully validated | CI import-check change (local only; never run on GitHub) |
| Planned, not implemented | splitting `mutation.py`, `cloning.py` and the other large modules; consolidating duplicated helpers; shorter dialogs; 85 % coverage |
| Blocked or deferred | remote CI (no remote); CRISPR/CFD and the AI Experiment Planner (owner decisions); licence; the two unreachable printers and legacy BLAST files (owner decisions) |

**Final verification.** Recorded in `docs/PHASE71_CONTINUATION.md`.
- Suite: 3.14.6 and 3.10.11.
- Golden master; Ruff; pyright gate.
- ReadOnly: 15/15 original files plus 27 internal-package files.
- `git status` clean after the final commit.

**Commits of the phase.**

| Commit | Content |
|---|---|
| `0bd2ace` | A: baseline, audit, characterization tests |
| `a45b8f2` | B: BLAST split |
| `c2dc998` | C1: explicit primer screen configuration |
| `61a8ba3` | C2: primer split |
| `d096e8e` | G: CI import check against the installed package |
| `f629add` | D: component tests |
| final | F/J: architecture documentation, reports |
