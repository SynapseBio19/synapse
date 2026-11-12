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
