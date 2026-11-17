# Phase 7 final report: engineering stabilization

Project Synapse, Phase 7, 2026-11-17. Working copy
`C:\Users\bayra\OneDrive\Desktop\Project Snyapse - Kopie`; the original
project was not modified. Governing document: the owner's Phase 7 Master
Execution Prompt (engineering quality only).

Binding restrictions, all kept:

- No new scientific features.
- CRISPR and CFD work excluded: `analyses/crispr.py`, `ui/crispr_menu.py`
  and their tests were not edited.
- The AI Experiment Planner stays **Coming Soon**.
- No licence chosen without the owner.

Status labels: **COMPLETED**, **PARTIALLY COMPLETED**, **DEFERRED** (owner
decision), **BLOCKED** (needs something outside this phase).

---

## A. Executive summary

Phase 7 made Synapse measurably more reliable to change. The scientific
behavior did not change; every refactoring was checked with byte-identical
output or the golden master.

| Measure | Start (`cf930cf`) | End | |
|---|---|---|---|
| Tests | 516 | 565 | all passing on Python 3.14.6 and 3.10.11 |
| Coverage (combined line/branch) | 75 % | 81 % | `main.py` 34 % -> 84 %, BLAST dialog 51 % -> 80 % |
| Pyright errors | 58 | 11 | the 11 are in the excluded CRISPR module; a CI gate forbids new ones |
| Ruff (default rules) | 20 findings, no configuration | 0, configured | one CRISPR-test finding excepted with the reason |
| FEATURE_STATUS rows contradicting code or tests | 13 known + 8 found | 0 found | now checked by a test |
| CI | none | GitHub Actions workflow, three jobs | has not run: no remote |
| `main()` length | 632 lines | 383 lines | output byte-identical |

Writing the new tests uncovered four defects, each now fixed and covered
by a test that failed before the fix:

1. **FEATURE_STATUS overstated testing.** It called eight functions
   "Execution-tested", but no test ran them.
2. **Remove Database left a file behind.** BLAST "Remove Database" left
   the version 5 metadata file.
3. **The BibTeX export broke LaTeX.** It produced entries on which LaTeX
   stops (`&`, `%` unescaped).
4. **Reports hid the genetic code.** The session report did not state
   the genetic code in its RNA and protein entries.

**Verdict.** Phase 7 is **PARTIALLY COMPLETED**:

- Documentation accuracy, type safety, CI configuration, dependency
  declaration, installation, error messages and report traceability are
  completed.
- Refactoring is partially completed: only `main()` was split. The large
  scientific modules were not split, deliberately (section E).
- The coverage target of 85 % was not reached (81 %).
- CI is configured but not operational, because there is no hosted
  repository.

---

## B. Initial engineering baseline

Measured 2026-11-17 at `cf930cf`, before any change. Recorded in
`docs/PHASE7_CONTINUATION.md`.

| Item | Value |
|---|---|
| Git | `main`, clean, 155 commits, no remote |
| Tests | 516, 0 failures, 5 skipped |
| Skip reasons | opt-in benchmark; four primer3-py comparisons, because Windows application control blocks its compiled module |
| Golden master | 12/12 |
| Pyright 1.1.414 at Python 3.10 | 58 errors |
| Pyright errors by location | tests 27; `crispr.py` 11; `literature.py` 8; `rna.py` 4; `orf.py` 3; `protein.py` 2; `core/report.py` 2; `blast.py` 1 |
| Ruff | pyflakes clean; 20 findings with the default rules; no configuration |
| Coverage | 75 % |
| Python | 3.14.6 (`.venv`), 3.10.11 (`.tools/`); no other interpreter |
| ReadOnly | 15/15 |
| Dependencies | `dependencies = []`; extras `reference`, `structure`, and `dev` (`dev` listed pytest, which is not used) |

**History sources.** `docs/PROJECT_HISTORY_AND_DECISIONS.md`, the Phase 6
reports and the two Claude Code sessions of Phase 6, read on 2026-11-13.
They have had no activity since. The claude.ai Project conversations and
the earlier "Claude Cowork" cycle are not accessible from this
environment and were not read.

---

## C. Documentation corrections (package A, `96d0f07`)

**Method.** Every one of the 271 FEATURE_STATUS rows was checked in two
ways:

- the functions it names were looked up in the code;
- their bodies were checked against line coverage of the full suite.

The first version of the check counted a function's `def` line, which
runs on import, as "executed". The corrected check counts only the body.

| Finding | Rows | Action |
|---|---|---|
| Out of date (named in the Phase 6 assessment) | 13 | Corrected after checking each against code and tests: RNA codon view (uses the ORF's genetic code since P1-9), S17, local blastp, the core API, custom BLAST databases, seven "Reachable" views that tests drive, the literature entry, the size guard |
| "Execution-tested" although no test ran the function | 8 | Tests written (`test_library_functions`, 14): double digest, `diagnose_thermodynamic_failure`, `normalize_indel_left`, `compare_proteins`, the Europe PMC, OpenAlex and Crossref providers |
| Named functions deleted in D16 part 2 | 2 | Names removed |
| Section headers with stale line counts | 14 | Counts removed |
| Validation ledger missing everything after 2026-11-11 | §17 | 24 validations listed; open list current |

**Prevention.** `tests/test_documentation_consistency.py` fails when the
inventory names code that does not exist.

**Earlier reports.** `PHASE6_FINAL_REPORT.md` and
`PHASE6_FINAL_ASSESSMENT.md` keep their text and get dated addenda; their
§12 status is to be read as PARTIALLY COMPLETE.

The other documents were updated with each package: README, METHODS,
AUDIT_STATUS, CHANGES, TERMINAL_UI_STANDARD and TEST_COVERAGE.

---

## D. Pyright and type safety (package B, `2f8f6e6`)

58 -> 11 errors. No check was disabled, no setting weakened and no broad
`Any` added.

| Root cause | Fix |
|---|---|
| A wrong annotation: `relevance_components: dict[str, float]` held lists and text | `TypedDict` with the real keys |
| Values guaranteed by the code but invisible to the checker: Shine-Dalgarno motif and position, mass-range bounds, BLAST runtime alias, date range, DOI | Narrowing that states the invariant (both bounds checked; `None` default) |
| Names bound on one branch only (ORF nesting) | One `nesting` tuple |
| An untyped helper (`get_rna_complement`) | `str -> str` |
| `hasattr` + attribute access in the report serializer | One `getattr` |
| Tests using results that can be `None` | `assert x is not None`; real `Message` headers; package import of a fixture; typed stand-ins |
| Biopython's stub declares `translate(table: str)` | One `pyright: ignore[reportArgumentType]` with the reason (the library documents and accepts int ids) |

**Remaining: 11 errors in `analyses/crispr.py`**, untouched by the owner's
decision. `scripts/ci/pyright_gate.py` allows exactly these and fails on
any other.

---

## E. Architectural refactoring (package C)

| Change | Commit | Evidence |
|---|---|---|
| `main()` 632 -> 383 lines: `compute_sequence_overview()` replaces 16 separate values passed twice through a 20-argument call; ORF step, protein step, no-ORF message and planner message are functions; the CRISPR branch is unchanged | `895fc5e` | Two scripted sessions byte-identical before and after (64,706 and 26,122 characters); golden master; session tests |
| Ruff default rule set clean (one-line statements, a redundant mid-module import, ambiguous names); `[tool.ruff]` configured | `c97695d` | Suite and golden master unchanged |

**Not changed, and why:**

- **The very large modules** (`primer.py` 5,183 lines, `mutation.py`
  4,670, `blast.py` 4,661, `cloning.py` 3,972) were not split. Splitting
  them is a multi-session change to scientific code. The task says not to
  split modules "arbitrarily" or to rewrite whole modules. The tests now
  give the coverage such a change needs. A plan follows below.
- **The self-tests** inside modules (`run_mutation_self_test`, 689 lines)
  stay where the menus call them.
- **Two unreachable report printers** (`print_protein_report`, about 600
  lines; `print_primer_report`) are kept. Deleting them is listed as an
  owner decision.

**Recommended next step.** Split the modules in this order, each in its
own commit, checked with the golden master and byte-identical rendered
menus:

1. move the self-tests to `analyses/selftests/` and re-export them;
2. split `blast.py` into client, local database manager, parsing and
   reports;
3. split `primer.py` into thermodynamics, screens, design and reporting.

The core/UI separation was verified again: `LayeringTests` pass, and core
imports and runs with `ui` unimportable.

---

## F. Testing and coverage (package D)

| New test file | Tests | What it covers |
|---|---|---|
| `test_library_functions` | 14 | the eight functions FEATURE_STATUS overclaimed |
| `test_documentation_consistency` | 2 | inventory names exist |
| `test_main_session` | 7 | `main.main()` in-process: file errors, the full DNA session, reports, every module entry except CRISPR, file and sequence switching, protein records, declared or cancelled types, invalid records |
| `test_blast_dialog` | 6 | result views after stubbed searches, exports, declined search, RID resume, database manager on a real `makeblastdb` database |
| `test_primer_and_overview_dialogs` | 4 | single-primer analysis, long-sequence overview |
| `test_literature_pipeline` | 8 | PubMed search, combined search with a failing provider, enrichment, five export formats |
| `test_message_rules` | 3 | no hand-padded fields, no prompts starting with a line break |
| `test_report_traceability` | 5 | provenance of every report entry; genetic code in ORF, RNA and protein |
| `test_blast_core_nt` (start of the day, Phase 6 close) | 4 | `core_nt` |

Totals: 516 -> 565 tests.

**Coverage by package** (coverage.py 7.10.7, combined line/branch):

| Package | Combined | Lines |
|---|---|---|
| `analyses` | 82.6 % | 86.6 % |
| `core` | 87.6 % | 90.7 % |
| `fileio` | 91.8 % | 94.8 % |
| `main.py` | 83.6 % | 84.3 % |
| `ui` | 74.2 % | 77.1 % |
| **Total** | **81 %** | |

**Why 85 % was not reached:**

- The CRISPR dialog (58 %, 200 lines missed) is excluded.
- `print_protein_report` (about 200 lines missed) is unreachable, and
  testing it would only raise the percentage.
- Remaining are paths that need a real NCBI database download or remote
  submission, and further restriction and cloning menu branches (71-75 %).

No test was weakened or removed.

---

## G. CI implementation (packages E/F, `37986d7`)

`.github/workflows/ci.yml` (GitHub Actions; no paid service) has three
jobs:

1. **`tests`.** Standard-library install on Python 3.10, 3.11, 3.12, 3.13
   and 3.14 (Ubuntu) and 3.14 (Windows). Steps: import check, golden
   master as its own step, full suite with skipped tests listed by reason.
2. **`reference`.** Biopython, primer3-py, ViennaRNA and BLAST+ installed.
   Only the opt-in benchmark may skip (`scripts/ci/check_skips.py
   --only`). A coverage report goes to the job summary.
3. **`lint`.** Ruff, and `scripts/ci/pyright_gate.py` in a `.venv` as
   `pyproject.toml` expects.

**Status: configured, not operational.** The repository has no remote,
so no run on GitHub exists. Local verification:

- the YAML was parsed and has 3 jobs, the matrix and the steps;
- every command was run;
- the reference gate fails here, as designed, because primer3-py is
  blocked on this machine.

---

## H. Dependencies and installation (package F, `37986d7`)

- **Runtime dependencies:** none (unchanged).
- **The `dev` extra** now lists the tools in use: pyright 1.1.414,
  Ruff 0.14.4 and coverage >= 7.10. pytest is gone: nothing imports it,
  no instruction used it, and the suite stays runnable with it. The
  unused pytest configuration is removed.
- **Python classifiers** for 3.10-3.14 are added.
- **Clean-environment check.** A new virtual environment (CPython 3.14.6)
  with `pip install .` only:
  - Synapse installs without dependencies;
  - its modules import from `site-packages`;
  - the golden master passes;
  - the suite passes with 22 visible skips.
- **Wheel check.** A wheel (`synapse-0.2.0.dev0-py3-none-any.whl`, 53
  files) was built and installed into a second clean environment; an ORF
  search ran from it, and ViennaRNA was reported as not available, which
  is correct.
- **Build artifacts** are ignored, and the README documents the
  development setup.
- **Python 3.11-3.13 were not run** here (no interpreters); they are in
  the CI matrix.

---

## I. Error handling and user experience (package G)

- **Escaping of the BibTeX export** (`ca5d53a`):
  - `&`, `%`, `$`, `#`, `_`, `~` and `^` are escaped in title, authors,
    journal, abstract and note;
  - the DOI is left as it is.
- **The BLAST database manager** removes the `.njs`/`.pjs` metadata file,
  and Database Information counts it (`cdbb7ac`).
- **Messages** (`3f48f1f`):
  - Five hand-padded single-primer lines use `style.field()`.
    TERMINAL_UI_STANDARD had wrongly said none remained; it is corrected.
  - Three BLAST prompts no longer start with a carriage return.
  - The BLAST export suffix error, the missing-cloning-module error and
    the literature date warning get their prefixes.
  - "No hits were returned" and the no-fallback remark are NOTEs.
  - `test_message_rules` keeps these forms from returning.
- **Five messages** in the mutation and literature dialogs were already
  fixed at the end of Phase 6 (`c24159f`).

---

## J. Reproducibility and traceability (package H, `22d0abb`)

A session report already recorded:

- in its header: software and Python version, platform, the input (hash,
  type and evidence) and its source file (path, format, hash);
- per entry: module, time, parameters, warnings, limitations and data
  hash.

New: the RNA entry states the genetic code and input mode, and the
protein entry the genetic code and reading frame. Both are read from the
results; nothing is inferred.

`test_report_traceability` runs a table-11 session and checks all of
this. METHODS "Reproducibility" describes what is recorded and where
external tools record their versions (BLAST, the enzyme database,
ViennaRNA).

---

## K. Remaining limitations and owner decisions

| Item | Status | Needed |
|---|---|---|
| Licence | BLOCKED (owner) | Choose a licence; confirm the copyright holder; optional `CITATION.cff` (`docs/LICENSING_AND_GOVERNANCE.md`) |
| Unused legacy files `blast_results/TEST_8_REFERENCE_blast.json`, `blast_test.json` | BLOCKED (owner) | Keep or delete |
| Unreachable `print_protein_report`, `print_primer_report` | BLOCKED (owner) | Keep or delete |
| CI on a hosted service | BLOCKED | A GitHub (or similar) remote |
| Splitting the large modules | PARTIALLY COMPLETED | Plan in section E |
| Coverage 85 % | PARTIALLY COMPLETED (81 %) | Further dialog tests (restriction, cloning, protein); CRISPR excluded |
| Pyright in `crispr.py` (11) | DEFERRED (owner: CRISPR excluded) | Lifting the exclusion |
| CFD scoring, other CRISPR items | DEFERRED (owner) | Lifting the exclusion |
| AI Experiment Planner | DEFERRED, **COMING SOON** | Owner specification |
| primer3-py comparisons on this machine | BLOCKED (Windows application control, not bypassed) | Run in CI or after the block is lifted |
| Python 3.11-3.13, Linux, macOS | Not run locally | The CI matrix, once it runs |
| Scientific limits from Phase 6 (KQ weight, genome-scale data, local swissprot, NEBcutter, experimental validation) | unchanged | See `docs/PHASE6_FINAL_REPORT.md` section K |

**AI Experiment Planner: COMING SOON — intentionally not implemented.**

---

## L. Final verification and Git state

Run on 2026-11-17 after the last code change:

| Check | Result |
|---|---|
| Full suite, CPython 3.14.6 | 565 run, 0 failures, 5 skipped (benchmark; four primer3-py comparisons blocked by application control) |
| Full suite, CPython 3.10.11 | 565 run, 0 failures, 23 skipped (Biopython 11, ViennaRNA 6, primer3-py 4, one 3.12-tokenizer test, benchmark) |
| Golden master | 12/12 |
| Pyright gate | 11 errors, all in the allow-list (`analyses/crispr.py`); none elsewhere |
| Ruff | all checks passed |
| Coverage | 81 % |
| Import check | `core.api`, `analyses.orf`, `analyses.rna_structure`, `ui.style`, from a clean install |
| Packaging | wheel built and installed in a clean environment; analysis ran from it |
| `git diff --check cf930cf..HEAD` | clean (CRLF files keep CRLF; `orf.py` keeps its original mixed endings) |
| ReadOnly | 15/15 restored |
| Background processes | none |

**Git.**

- `main`, working tree clean after the commit of this report.
- Phase 7 commits: `96d0f07`, `2f8f6e6`, `cdbb7ac`, `bedcd42`, `ca5d53a`,
  `c97695d`, `895fc5e`, `37986d7`, `3f48f1f`, `22d0abb`, `236008f`, and
  the commit of this report.
- No history rewrite, reset or force-push.

**Assessment.** Synapse is a cleaner, better-tested and better-documented
research and teaching platform than at the start of the phase:

- its documentation is now checked against its code;
- type errors outside the excluded module are gone;
- its session code is shorter;
- it can be installed and checked reproducibly.

It is not yet professionally engineered in every respect. Large modules
remain, CI has not run on a server, no licence exists, and the scientific
limits recorded in Phase 6 are unchanged.
