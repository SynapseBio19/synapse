# Phase 7.1 architecture audit

Written on 2026-11-18, before any Phase 7.1 refactoring, at HEAD `af7cd0a`.
The numbers come from an AST scan of `main.py`, `analyses/`, `core/`,
`ui/`, `fileio/` and `data/` (CRISPR code is measured, not examined), from
the import graph, and from reading the BLAST and primer modules.

## 1. Verified baseline

| Item | Value (2026-11-18) |
|---|---|
| Branch, HEAD | `main`, `af7cd0a`, 167 commits, working tree clean, no Git remote |
| Python | CPython 3.14.6 (`.venv`), CPython 3.10.11 (`.tools/python-3.10.11`) |
| Tests, 3.14.6 | 565 run, 0 failures, 1 skipped (the opt-in benchmark) |
| Tests, 3.10.11 | 565 run, 0 failures, 23 skipped (optional packages not installed there) |
| primer3-py | **loads on this machine today.** The four primer3-py comparisons ran and passed (`tests.test_primer_thermodynamics`, `tests.test_primer_structure_reference`). On 2026-11-17 Windows application control blocked its compiled module and they skipped. |
| Golden master | 12/12 |
| Ruff | clean |
| Pyright | 11 errors, all in `analyses/crispr.py` (on the allow-list of `scripts/ci/pyright_gate.py`) |
| Coverage | 81 % combined line/branch (18,008 statements, 2,787 missed); `analyses/blast.py` 77 %, `analyses/primer.py` 86 % |
| ReadOnly | 15/15 protected files set |

The Phase 7 report's values are confirmed, except that the skip count on
3.14.6 dropped from 5 to 1 because primer3-py now loads.

## 2. Size and shape

Module sizes are physical lines. "Self-test" is the number of lines in
built-in `run_*self_test` and `_test_*` functions.

| Module | Lines | Functions | Classes | Self-test lines |
|---|---|---|---|---|
| `analyses/primer.py` | 5,182 | 103 | 16 | 579 |
| `analyses/mutation.py` | 4,670 | 73 | 15 | 706 |
| `analyses/blast.py` | 4,667 | 150 | 15 | 654 (plus a 54-line sample XML) |
| `analyses/cloning.py` | 3,972 | 57 | 0 | 239 |
| `analyses/protein.py` | 2,544 | 32 | 0 | 0 |
| `analyses/variant_analysis.py` | 2,531 | 72 | 12 | 261 |
| `ui/protein_menu.py` | 1,833 | 24 | 0 | 0 |

Largest functions (lines):

| Function | Lines | Note |
|---|---|---|
| `analyses/mutation.py: run_mutation_self_test` | 689 | self-test |
| `ui/protein_menu.py: print_protein_report` | 599 | unreachable (owner decision) |
| `analyses/orf.py: find_orfs` | 593 | scientific core |
| `ui/blast_menu.py: blast_menu` | 571 | dialog |
| `ui/rna_menu.py: print_rna_report` | 566 | report |
| `ui/report.py: print_results` | 553 | report |
| `analyses/cloning.py: evaluate_cloning_strategy` | 425 | scientific core |
| `analyses/primer.py: run_thermodynamic_self_test` | 379 | self-test |
| `analyses/blast.py: LocalBlastEngine.run` | 247 | local BLAST execution |
| `analyses/blast.py: LocalBlastDatabaseManager.scan_database` | 197 | database discovery |

## 3. Import graph

- **No circular imports at module level.** Module-level imports point only
  from `ui` to `analyses`/`core`/`data`, from `analyses` to `core`/`data`,
  and from `core` to `core`/`data`. `LayeringTests` already enforce that
  `analyses` and `core` never import `ui` (the `ui` imports in `analyses`
  sit inside `if __name__ == "__main__":` blocks).
- `core/api.py` and `core/pipeline.py` import the analyses lazily inside
  functions. This keeps the start-up light and is intended.
- `analyses/blast.py` imports `fileio.loader` and `core.orf_context` lazily
  inside two functions.

## 4. Findings

Priority: **P1** high architectural impact and handled in Phase 7.1,
**P2** worth doing next, **P3** recorded only. "Risk" is the regression
risk of the proposed change.

### F1 (P1). `analyses/blast.py` holds ten unrelated responsibilities

1. **Component.** The whole module, 4,667 lines.
2. **Current responsibility.** By line range:
   - configuration and catalogs: endpoints, rate limits, programs, databases, local catalog (1-360);
   - data model, 10 dataclasses (360-690);
   - query validation and FASTA loading (690-930);
   - NCBI URL API client (933-1247, 315 lines);
   - local database model, executable discovery and the database manager (1248-2514; `LocalBlastDatabaseManager` alone is 985 lines);
   - runtime alias views and the BLAST environment (2515-2750);
   - local BLAST execution (2751-3012);
   - BLAST XML parsing (3013-3283);
   - evidence scoring, ranking and alignment text (3284-3462);
   - the ORF adapter and the `run_blast_analysis` workflow (3463-3652);
   - serialization, exports, text report and run comparison (3653-3873);
   - the built-in self-test (3874-4667, 794 lines).
3. **Problem.** Every change, whether to the XML parser, the HTTP client or
   the database manager, touches the same file. The parser and the
   evidence score can only be read together with HTTP, `subprocess` and
   file-system code. Reviewing a parser change means scrolling past 2,000
   lines of database management. A sixth of the module is test code that
   ships with the application.
4. **Proposed improvement.** Move the implementation into an internal
   package `analyses/_blast/` with one module per responsibility.
   `analyses/blast.py` stays as the public module: it is the protected
   file, the place every caller imports from, and the `python -m` entry
   point. It re-exports every name it had.
   - A package named `analyses/blast/` cannot sit next to `blast.py`, and
     replacing the protected file would break the protection rule and
     every `from analyses import blast`.
5. **Risk.** Medium.
   - `BlastRun.ranked_hits` calls `rank_hits`, and the evidence code needs
     `BlastRun`. If ranking lands in a different module from the model,
     the two modules import each other.
   - Three functions locate files through `Path(__file__)` (database root,
     executable search, user paths). In a subpackage they would silently
     point one level too deep.
   - Two tests patch module attributes (`B.urlopen`, `B.time.sleep`). After
     the move these must name the module that calls them, or the patch no
     longer takes effect.
6. **Dependencies affected.**
   - `ui/blast_menu.py` imports 36 names, two of them private
     (`_evalue_text`, `_infer_scientific_name_from_definition`).
   - `core/api.py` uses two names.
   - 11 test modules and the golden-master scenario.
   - `tests/test_builtin_selftests.py` calls `analyses.blast.run_self_test`.
   - CRISPR does not import BLAST.
7. **Tests required.**
   - Namespace test: all 125 top-level names stay reachable.
   - Characterization snapshot: results, report text and the CSV/JSON/text
     exports; database root and executable search paths.
   - The existing BLAST tests, rendered menus, golden master and session
     tests.
8. **Expected benefit.** Each responsibility can be read, reviewed and
   tested on its own. The parser and the evidence score no longer pull in
   network or process code. The test code is separated from the
   application code.

### F2 (P1). `analyses/primer.py` is a stack of historical layers

1. **Component.** The whole module, 5,182 lines.
2. **Current responsibility.** By line range:
   - sequence helpers and the legacy design engine with its data model (1-1075);
   - a self-test (1075-1180) and an `if __name__` block;
   - the v2/v3 nearest-neighbor thermodynamics (1256-2280);
   - the professional data model and scoring (2281-3008);
   - two self-tests (3009-3424);
   - preset conditions and the public API list (3425-3586);
   - a second `if __name__` block;
   - the backward-compatible `calculate_tm_nearest_neighbor` with four overloads (3587-3826);
   - "v4"-"v11" sections: the candidate screen, structure checks, pair
     ranking, the professional entry point and restriction-site design
     (3861-5182).
3. **Problem.** Reading order does not match dependency order.
   - The legacy engine at line 566 calls `calculate_tm_nearest_neighbor`,
     defined 3,100 lines later. `analyze_hairpin` (line 4393) is used from
     line 602 on.
   - This works only because Python resolves module globals at call time.
   - `GAS_CONSTANT` is assigned twice, with the same value.
   - The "v5"-"v11" prefixes describe history, not alternatives. All of
     them are live code.
   - Self-tests are spread over four places.
4. **Proposed improvement.** An internal package `analyses/_primer/`.
   Modules follow the dependency order: sequence helpers, model,
   thermodynamics, structure, specificity, legacy design, professional
   model and scoring, candidate screen, pair ranking, professional
   workflow, self-tests. `analyses/primer.py` stays as the public module,
   with both `__main__` blocks.
5. **Risk.** Medium to high.
   - Forward references must become imports without a cycle. Checked: the
     proposed order is acyclic.
   - Three `lru_cache` wrappers must stay single objects.
   - The hidden global (F3) has to move into one module.
6. **Dependencies affected.**
   - `ui/primer_menu.py`, `core/api.py` and `core/pipeline.py`.
   - 11 test modules besides the golden master, among them the
     screen-equivalence tests and the primer3-py comparisons.
7. **Tests required.**
   - Namespace test for 140 names.
   - Characterization snapshot:
     - legacy design and its summaries;
     - professional design under all three Tm models and high salt;
     - restriction-site design;
     - single-primer analyses;
     - Tm comparisons and structure checks for four primers;
     - three condition sets;
     - all four self-tests.
   - The existing primer suites and the primer3-py comparisons (run today).
8. **Expected benefit.** Each module can be read top-down without forward
   references. Thermodynamics can be changed and reviewed without the
   design engine. The candidate screen and the ranking become separate,
   testable units.

### F3 (P1). Hidden global configuration in the primer candidate screen

1. **Component.** `_V5_ACTIVE_CONFIG` in `analyses/primer.py`.
2. **Current responsibility.**
   - `_v7_generate_professional_primer_candidates` writes the
     configuration of the current design into this global.
   - `_v5_candidate_quality` reads the optimum Tm and GC from it.
3. **Problem.**
   - A result depends on the last call, not only on the arguments.
   - Two tests set the global from outside before calling internal
     functions.
   - The state is not thread-safe.
   - A module split would split the variable: a test assigning
     `primer._V5_ACTIVE_CONFIG` would no longer reach the screen.
4. **Proposed improvement.** Pass the configuration as an argument. Both
   call sites (`_v7_screening_score`, `_v7_measure_candidates`) already
   receive it as their `config` parameter. In every public path it is the
   same object the global holds.
5. **Risk.** Low. The characterization snapshot covers all three Tm models,
   and the screen-equivalence tests compare the screen with the exhaustive
   measurement.
6. **Dependencies affected.** Two tests (`test_primer_corrections`,
   `test_primer_screen_equivalence`) that set the global.
7. **Tests required.** Snapshot, screen-equivalence tests, and a new test
   that the score depends on the configuration that is passed.
8. **Expected benefit.** The screen becomes a pure function of its
   arguments.

### F4 (P2). Self-tests inside production modules

- **Where.** `mutation.py` (706 lines), `blast.py` (654 + sample XML),
  `primer.py` (579), `literature.py` (269), `variant_analysis.py` (261),
  `cloning.py` (239), `restriction.py` (202), `sequence.py` (136).
- **Problem.** They are part of the menus' "self-test" entries, so they
  must stay callable. But they double the reading load of the modules.
- **Proposal.** Move them next to the implementation. Done for BLAST and
  primer in Phase 7.1 as part of F1/F2; the others follow the same
  pattern.
- **Risk.** Low. `tests/test_builtin_selftests.py` runs every one.

### F5 (P2). `mutation.py` (4,670) and `cloning.py` (3,972) not examined in depth

- They are the next-largest modules.
- `cloning.py` has no classes. `evaluate_cloning_strategy` (425 lines) and
  `automatic_cloning_planner` (208) carry most of the logic, and their
  coverage is 71-75 % in the menus.
- Phase 7.1 does not touch them. The BLAST and primer pattern (public
  module plus internal package, a namespace test and a snapshot) applies
  to them unchanged.

### F6 (P2). Duplicated sequence helpers

- **Where.**
  - `reverse_complement` is defined seven times (`core/seqtools.py`,
    `core/orf_context.py`, cloning, crispr, mutation, primer, restriction).
  - `normalize_sequence` four times.
  - FASTA is read by `fileio.loader`, `blast.load_fasta_queries` and
    `crispr.load_fasta`.
- **Problem.** The copies differ in validation and error messages (IUPAC
  handling, `ValueError` versus module-specific errors). A fix in one does
  not reach the others.
- **Proposal.** Consolidate on `core.seqtools` one function at a time.
  Each step first records the old error messages and IUPAC behavior.
- **Risk.** Medium: error texts are part of the rendered menus.
- Not done in Phase 7.1, because the gain is smaller than the regression
  risk under a "no behavior change" rule.

### F7 (P2). Very long dialog and report functions

- **Where.** `blast_menu` (571 lines), `print_rna_report` (566),
  `print_results` (553), `print_orf_report` (405),
  `run_mutation_analysis_ui` (395), `main()` (383); `crispr_menu` (496) is
  excluded.
- **Problem.** These are long `if`/`elif` chains over menu entries.
- **Proposal.** Extract one function per entry, as was done for `main()` in
  Phase 7. The rendered-menu tests compare byte for byte.
- Not done here, because the UI is not to be redesigned and the gain is
  local.

### F8 (P3). Unreachable printers

`print_protein_report` (599 lines) and `print_primer_report` have no
caller. They stay until the owner decides
(`docs/LICENSING_AND_GOVERNANCE.md` §4).

### F9 (P3). Configuration resolved at import time

- **What.** `analyses.blast` reads `SYNAPSE_NCBI_TOOL`, `SYNAPSE_NCBI_EMAIL`
  and `SYNAPSE_BLAST_DB_ROOT` when it is imported. `DEFAULT_LOCAL_DB_ROOT`
  also inspects the file system at that point.
- **Problem.** A change to the environment after import has no effect.
- **Decision.** Kept. Making it lazy would change when the database root
  is decided, which is a behavior change. The split moves this code
  unchanged.

### F10 (P3). Tests coupled to module layout

- **What.**
  - `mock.patch.object(B, "urlopen")` and `mock.patch.object(B.time, "sleep")`
    patch the module that happens to contain the client.
  - The two primer tests from F3 write a module global.
- **Problem.** Such tests break or, worse, silently stop patching when code
  moves.
- **Proposal.** Patch the module that calls the function.
  `test_module_namespaces` also guards the names.

### F11 (P3). Private helpers used across layers

- **What.** `ui/blast_menu.py` uses `_evalue_text` and
  `_infer_scientific_name_from_definition`.
- **Decision.** Kept as they are. Renaming them would be an API change. The
  public module re-exports them.

## 5. Plan and order

1. Characterization before refactoring:
   - `tests/test_module_namespaces.py` (264 names);
   - `tests/test_characterization.py` with the snapshot
     `tests/fixtures/characterization/blast_primer.json` (78 scenarios,
     written from `af7cd0a`).
2. F1, BLAST: create `analyses/_blast/`. The code moves unchanged, except:
   - `Path(__file__)` lookups go through a constant for the `analyses`
     directory, so every path stays the same;
   - the two tests from F10 patch the client module.
3. F3, then F2, primer: remove the hidden global first, in its own commit
   with its own test. Then move the code.
4. Documentation: `docs/PHASE71_ARCHITECTURE.md`, then the final report.

Moves are done by a script that cuts the module at statement boundaries
and generates the imports from the names each piece uses. Moved code is
therefore byte-identical apart from the documented path change. The
script's dependency check confirms that the new packages have no import
cycle.
