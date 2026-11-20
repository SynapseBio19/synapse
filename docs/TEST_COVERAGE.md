# Test coverage

Measured on 2026-11-12 at the commit that adds this file (Phase 6, task
T7). Coverage shows which code the automated suite executes. It does not
show that the executed code is correct; that is the job of the assertions
and of the external reference comparisons (`docs/METHODS.md`).

## Method

- Tool: coverage.py 7.10.7 (C extension), installed into the project's
  `.venv` on 2026-11-12 under the approval policy in
  `docs/PROJECT_SYNAPSE_COMPLETE_HANDOVER.md`, Part X. It is a development
  tool, not a runtime dependency.
- Interpreter: CPython 3.14.6, Windows 11.
- Command (from the project root):

  ```bash
  .venv/Scripts/python.exe -m coverage run --branch --source=analyses,core,ui,fileio,data,main -m unittest discover -s tests -t .
  .venv/Scripts/python.exe -m coverage report --sort=cover
  ```

- Line and branch coverage; "Total" is coverage.py's combined figure.

## Result

Suite: 450 tests, all passing (1 opt-in benchmark skipped).

| Package | Statements | Lines covered | Branches | Branches covered |
|---|---:|---:|---:|---:|
| `analyses` | 11,547 | 82 % | 4,360 | 66 % |
| `core` | 1,114 | 90 % | 300 | 76 % |
| `data` | 21 | 76 % | 0 | – |
| `fileio` | 213 | 95 % | 106 | 86 % |
| `ui` | 4,481 | 53 % | 1,532 | 44 % |
| `main.py` | 299 | 34 % | 86 | 23 % |
| **Total** | 17,675 | 74 % | 6,384 | 61 % (combined 70.9 %) |

The first measurement on the same day, before the tests below were added,
gave 68 % combined.

Per file (combined figure): `core/alphabets.py` 99 %, `core/api.py` 93 %,
`core/report.py` 91 %, `analyses/crispr.py` 89 %, `analyses/orf.py` 88 %,
`analyses/sequence.py` 85 %, `analyses/primer.py` 84 %,
`analyses/protein.py` 84 %, `analyses/mutation.py` 80 %,
`analyses/variant_analysis.py` 79 %, `analyses/cloning.py` 75 %,
`analyses/rna.py` 75 %, `analyses/blast.py` 74 %, `core/pipeline.py` 74 %,
`analyses/restriction.py` 70 %, `analyses/literature.py` 53 %.

## What the first measurement found

- The protein and RNA dialogs ran only to their first menu (9 % and
  14 %); no view behind the menus was ever drawn. `tests/test_execution_coverage.py`
  now opens every view and checks its headings and content (protein
  47 %, RNA 75 % afterwards).
- Eleven science functions were never called, among them the BLAST JSON,
  CSV and text exports, the CRISPR CSV and FASTA exports, the base- and
  prime-editing reports, and the variant grouping, density, formatting
  and report helpers. The same test file runs each with checks of the
  output.
- Running the BLAST text export showed a defect (D24): the report ends its
  lines with CR LF and was written in text mode, so every line ended in
  CR CR LF. The local BLAST query file had the same fault. Both are now
  written with `newline=""`.

## Known limits of the figures

- `main.py` and `ui/report.py` are mostly executed by `tests/test_main_smoke.py`,
  which runs `main.py` in a subprocess with scripted input. coverage.py does
  not follow subprocesses in this configuration, so their figures are lower
  bounds; the smoke test does draw the sequence overview, the ORF, RNA and
  protein screens and the pipeline status.
- Remaining uncovered code is mostly behind network or tool access that the
  suite deliberately does not use: live NCBI BLAST and the literature
  providers (`analyses/literature.py` 53 %, `ui/literature_menu.py` 23 %),
  BLAST database downloads and `makeblastdb`, local BLAST execution, and
  file-export dialogs (`ui/variant_analysis_menu.py` 21 %,
  `ui/mutation_menu.py` 38 %, `ui/orf_menu.py` 35 %).
- 100 % coverage is not a goal in itself; the gaps above are recorded so
  that a reader knows which dialogs have only been rendered, not driven to
  the end.

## Remeasured 2026-11-17 (continuation run)

Same tool and command. Suite: 516 tests, all passing (5 skipped: the
opt-in benchmark and four primer3-py comparisons, because Windows
application control blocked primer3-py's compiled module that day).

`tests/test_dialog_paths.py` (17 tests) drives the four dialogs that were
least covered to the end of their branches: operations before any data
exist, missing files, empty and invalid answers, cancel, settings,
selection, and the VCF, TSV and literature exports into a temporary
folder (the literature providers are stubbed; nothing is sent).

| File | 2026-11-12 | 2026-11-17 |
|---|---:|---:|
| `ui/variant_analysis_menu.py` | 37 %* | 87 % |
| `ui/orf_menu.py` | 35 % | 76 % |
| `ui/mutation_menu.py` | 38 % | 72 % |
| `ui/literature_menu.py` | 23 % | 72 % |
| **Total (combined)** | 72 %* | 75 % |

\* measured on 2026-11-17 before the new tests (the 2026-11-12 table above
gives 21 % for the variant dialog and 70.9 % in total; tests added since
raised them).

Found while writing them, and corrected: three literature messages and
two mutation messages were printed without their ERROR/NOTE/WARNING
prefix, and the ORF settings note did not say that ORFfinder's minimal
length excludes the stop codon.

Still mostly uncovered: `ui/report.py` (12 %, driven only through the
`main.py` subprocess, which coverage does not follow), the protein and
BLAST dialogs (47 % and 51 %), and the local BLAST database download and
`makeblastdb` paths. The CRISPR dialog was deliberately not extended in
this run.


## Remeasured 2026-11-18 (end of Phase 7.1)

Same tool and command. Suite: 601 tests, all passing (1 skipped: the
opt-in benchmark; primer3-py loads again, so its comparisons ran).

Since Phase 7.1 `analyses/blast.py` and `analyses/primer.py` keep their
code in `analyses/_blast/` and `analyses/_primer/`, so coverage is
reported per responsibility. `tests/test_blast_components.py` and
`tests/test_primer_components.py` were written for the parts that the
split showed as least reached.

| Scope | 2026-11-17 | 2026-11-18 |
|---|---:|---:|
| BLAST (`blast.py`, then `blast.py` + `_blast/`) | 77 % | 84.9 % |
| Primer (`primer.py`, then `primer.py` + `_primer/`) | 86 % | 91.7 % |
| `analyses` | 82.6 % | 84.6 % |
| `core` | 87.6 % | 87.6 % |
| `ui` | 74.2 % | 74.2 % |
| `fileio` | 91.8 % | 91.8 % |
| `main.py` | 83.6 % | 83.6 % |
| **Total** | **81 %** | **82.3 %** |

Least-covered BLAST modules now:
- `_blast/runtime` (67 %): the cleanup after a failed alias view;
- `_blast/ncbi` (72 %): HTTP retry branches;
- `_blast/workflow` (73 %);
- `_blast/databases` (75 %): the NCBI download and verification against a
  real database.

The 85 % target is not reached: the CRISPR dialog stays excluded, and
`print_protein_report` stays unreachable until the owner decides.

## Remeasured 2026-11-18 (Phase 8)

Same tool and command. Suite: 608 tests, all passing. 7 skipped: the
opt-in benchmark and six ViennaRNA tests, because Windows application
control blocked ViennaRNA's DLL that day.

| Scope | Phase 7.1 end | Phase 8 |
|---|---:|---:|
| Mutation (`mutation.py` + `_mutation/`) | – | 88.8 % |
| Cloning (`cloning.py` + `_cloning/`) | – | 85.0 % |
| `analyses` | 84.6 % | 85.5 % |
| `core` | 87.6 % | 87.0 % |
| `ui` | 74.2 % | 77.1 % |
| `fileio` | 91.8 % | 91.8 % |
| `main.py` | 83.6 % | 83.6 % |
| **Total** | **82.3 %** | **83.7 %** |

- `ui` rose because the 638 unreachable lines of the two deleted report
  printers no longer count, not because of new tests.
- `core` fell by 0.6 points. The likely cause is the ViennaRNA entry
  point in `core/api.py`, which only the skipped tests reach; this was not
  checked line by line.

## Remeasured 2026-11-19 (Phase 9)

Same tool and command, with `synapse_gen` added to the measured packages.
Suite: 659 tests, all passing. 5 skipped: the opt-in benchmark and the
four opt-in refseq_rna searches. The primer3-py and ViennaRNA comparisons
ran.

| Scope | Phase 8 | Phase 9 |
|---|---:|---:|
| `analyses` | 85.5 % | 87.3 % |
| `core` | 87.0 % | 87.7 % |
| `ui` | 77.1 % | 77.7 % |
| `fileio` | 91.8 % | 91.8 % |
| `synapse_gen` (the session, formerly `main.py` at 83.6 %) | – | 83.0 % |
| **Total** | **83.7 %** | **85.0 %** |

The 85 % target of Phase 7 is reached. It is reached by tests of real
behavior, not by tests written for the percentage:
- component tests;
- the genetic-code, ORF-differential and determinism tests;
- the download check;
- the command and path tests.

## Remeasured 2026-11-20 (Phase 10)

Same tool and command (`--branch`, sources `analyses,core,ui,fileio,data,synapse_gen`),
CPython 3.14.6. Suite: 722 tests, all passing, 5 opt-in skips. For a fair
comparison the code of `a2898aa` (end of Phase 9) was measured again with
the same command on the same day.

| Scope | `a2898aa` | Phase 10 |
|---|---:|---:|
| `ui/cloning_menu.py` | 74 % | 82 % |
| `ui/restriction_menu.py` | 75 % | 85 % |
| `ui/protein_menu.py` | 83 % | 87 % |
| `analyses/_literature/` (new) | – | 81-100 % per module |
| `core/windows.py` (new) | – | 100 % |
| **Total (combined)** | **85.0 %** | **86.1 %** |

The gain comes from the dialog scenarios of
`tests/test_protein_restriction_cloning_dialogs.py` and the differential,
benchmark and user-guide tests; no test was written for the percentage.
