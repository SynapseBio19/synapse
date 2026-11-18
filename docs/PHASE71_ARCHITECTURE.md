# Synapse architecture after Phase 7.1

This describes the code as it is at the end of Phase 7.1 (2026-11-18). The
findings that led to the changes are in `docs/PHASE71_ARCHITECTURE_AUDIT.md`;
the change records are in `docs/CHANGES.md` (section "Phase 7.1").

## 1. Layers

```
main.py                 terminal session: sequence loading, the module menu
ui/                     dialogs and reports (draw through ui/style.py)
core/                   data model, pipeline and provenance, public API, alphabets
analyses/               the scientific modules (one public module per analysis)
  _blast/, _primer/     internal implementation packages (Phase 7.1)
fileio/                 FASTA and GenBank reading
data/                   genetic codes, enzymes, amino-acid data
```

- **Import rules.** Enforced by `tests/test_terminal_ui.py`, LayeringTests.
  - `analyses` and `core` never import `ui`, except inside
    `if __name__ == "__main__":`.
  - `ui` uses `analyses`, `core` and `data`. `core/api.py` and
    `core/pipeline.py` import the analyses lazily inside functions.
  - At module level there is no import cycle.
- **Since Phase 7.1.**
  - The scan behind these rules (and the spelling, vocabulary, notation and
    feature-inventory checks) covers `analyses/**`, so the internal
    packages are held to the same rules.
  - No internal module imports a public one; `test_module_namespaces`
    checks this.

## 2. Public modules and internal packages

`analyses/blast.py` and `analyses/primer.py` remain the only modules that
callers import:
- `ui/blast_menu.py`, `ui/primer_menu.py`, `core/api.py` and
  `core/pipeline.py`;
- the tests;
- `python -m analyses.blast` and `python -m analyses.primer`.

Each public module:

- keeps its original docstring, which describes the scientific scope and
  limits;
- re-exports every name it defined at `9077f1e`, private helpers included
  (`ui/blast_menu.py` uses two of them). The only exception is the removed
  global `_V5_ACTIVE_CONFIG`;
- keeps its `if __name__ == "__main__":` block or blocks.

The implementation lives in `analyses/_blast/` and `analyses/_primer/`.
The leading underscore marks them as internal: their module layout may
change without notice, the public module's names may not.

- **Why two directories.** A package `analyses/blast/` cannot coexist with
  `analyses/blast.py`. Replacing the module would have removed a file
  under ReadOnly protection and changed every import.
- **Patching in tests.** To patch something a component calls (for
  example `urlopen`, `time.sleep`, `subprocess.run`), patch it in the
  module that calls it: `analyses._blast.ncbi`,
  `analyses._blast.local_engine`. Patching the public module has no
  effect on the implementation.

## 3. BLAST (`analyses/_blast/`)

| Module | Lines | Responsibility | Depends on |
|---|---|---|---|
| `paths` | 34 | `ANALYSES_DIR` (the `analyses` directory), resolution of user-supplied paths | – |
| `config` | 244 | NCBI endpoint, rate limits and timeouts (from `SYNAPSE_NCBI_*` at import), programs, remote databases, the local database catalog, the default database root (from `SYNAPSE_BLAST_DB_ROOT` or the project's `data/blast_db`) | `paths` |
| `models` | 428 | dataclasses (`BlastQuery`, `BlastParameters`, `BlastHSP`, `BlastHit`, `BlastRun`, `BlastRunMetadata`, `HitEvidence`, `BlastComparison`, `LocalBlastDatabase`, `SelfTestResult`), `BlastAnalysisError`, `rank_hits` | `config` |
| `query` | 296 | sequence normalization and type inference, program/query compatibility, FASTA queries, the ORF adapter | `config`, `models`, `paths` |
| `parsing` | 323 | BLAST XML (outfmt 5) into the model; taxonomy fallback from the hit definition | `config`, `models` |
| `ncbi` | 361 | `NCBIBlastClient`: submit, poll within the NCBI limits, fetch, resume an RID; transient-error retries | `config`, `models`, `parsing`, `query` |
| `executables` | 145 | `BlastExecutableManager`: finding BLAST+ and `update_blastdb.pl` | `config`, `models`, `paths` |
| `runtime` | 232 | command references, temporary alias views for nested volumes, the `BLASTDB` environment | `models` |
| `databases` | 1,121 | `LocalBlastDatabaseManager`: discovery, verification, manifest, NCBI download, custom databases, removal | `config`, `executables`, `models`, `runtime` |
| `local_engine` | 302 | `LocalBlastEngine`: one local search, the short-query rescue | `config`, `executables`, `models`, `parsing`, `query`, `runtime` |
| `evidence` | 166 | evidence score and profile, taxonomy grouping, alignment text | `config`, `models` |
| `reporting` | 232 | dictionaries, JSON/CSV/text exports, the text report, run comparison | `evidence`, `models`, `paths` |
| `workflow` | 118 | `run_blast_analysis`: one search, remote or local, from a sequence or an ORF | `config`, `databases`, `local_engine`, `models`, `ncbi`, `query` |
| `selftest` | 841 | the built-in self-test and its sample XML | all of the above |

**Data flow of a search.**

1. `query.build_query` / `build_orf_query` -> `BlastQuery`.
2. Remote search: `ncbi.NCBIBlastClient.run`, or local search:
   `databases.LocalBlastDatabaseManager` + `local_engine.LocalBlastEngine.run`.
3. Both produce outfmt-5 XML, which `parsing.parse_blast_xml` turns into
   a `BlastRun`.
4. `evidence` scores and ranks the hits; `reporting` exports them.
5. `ui/blast_menu.py` draws every step. `workflow.run_blast_analysis` is
   the same chain without a dialog (`core/api.py`).

**Design decisions.**
- `rank_hits` sits in `models` because `BlastRun.ranked_hits` calls it.
- `LocalBlastDatabase` is in `models` because `runtime` and `databases`
  both need it.
- `LocalBlastDatabaseManager` stays one class of 985 lines. Its methods
  share the manifest and the volume layout, and splitting a class would
  have changed its interface.

## 4. Primer (`analyses/_primer/`)

| Module | Lines | Responsibility | Depends on |
|---|---|---|---|
| `sequence` | 116 | DNA normalization (cached), target validation, reverse complement, GC, clamp, repeats | – |
| `models` | 146 | legacy configuration and results (`PrimerDesignConfig`, `PrimerCandidate`, `PrimerPair`, `PrimerDesignReport`, ...) | `sequence` |
| `thermodynamics` | 1,581 | Wallace and salt-adjusted Tm; SantaLucia 1998/2004 and the legacy parameter set; salt and concentration terms; duplex thermodynamics (cached); `calculate_tm_nearest_neighbor` and its overloads; condition presets and summaries | `models`, `sequence` |
| `structure` | 194 | hairpin, self-dimer and heterodimer screens (sequence-level heuristics) | `models`, `sequence` |
| `specificity` | 91 | screen against caller-supplied references, binding sites of a pair | `models`, `sequence` |
| `design` | 589 | the legacy engine: candidate scoring and generation, pairs, `design_primers`, `analyze_primer` | `models`, `sequence`, `specificity`, `structure`, `thermodynamics` |
| `professional` | 685 | `ProfessionalPrimerConfig`, advanced candidates and pairs, thermodynamic scores, `analyze_primer_professional`, failure diagnosis | `models`, `sequence`, `structure`, `thermodynamics` |
| `screening` | 425 | the professional candidate pool: all windows measured without hard filters, the screening score, the Tm-tier screen with its upper bound | `design`, `models`, `professional`, `sequence`, `thermodynamics` |
| `ranking` | 447 | pair quality, diversity selection, ranking, pair summary | `design`, `professional` |
| `workflow` | 277 | `design_primers_professional`, restriction-site design | `professional`, `ranking`, `screening`, `sequence` |
| `selftest` | 745 | the four built-in self-tests, `PRIMER_V3_PUBLIC_API` | all of the above |

**Data flow of a professional design** (`ui/primer_menu.py` and
`core/api.py`):

1. `workflow.design_primers_professional`
2. -> `screening.generate_professional_primer_candidates`, per
   orientation. It measures each window with `design._make_candidate`
   (Tm from `thermodynamics`, structures from `structure`) and scores it
   with `_screening_quality(item, config)`.
3. -> `ranking.rank_professional_primer_pairs`: pairs from
   `design.build_primer_pair`, with the specificity screen when references
   are given.
4. -> a `ProfessionalPrimerDesignReport` with failure reasons and
   warnings.

**Design decisions.**
- The screening score takes the configuration as an argument (Phase 7.1,
  F3). There is no module-level state left apart from three `lru_cache`
  wrappers. Those are pure caches of normalization, nearest-neighbor sums
  and duplex thermodynamics.
- The `_v5`-`_v11` prefixes were kept. Renaming would change 21 private
  names that tests and the namespace contract use. The package docstring
  and the section comments say that they are history, not versions.
- The legacy engine (`design`) and the professional engine coexist, as
  before. Both are public API.

## 5. Testing strategy for refactorings

| Layer | Test | Guards |
|---|---|---|
| Names | `tests/test_module_namespaces.py` | every name of `analyses.blast` (125) and `analyses.primer` (139) stays reachable; each internal module imports alone and complete; public names are the implementation objects |
| Results | `tests/test_characterization.py` + `tests/fixtures/characterization/blast_primer.json` | 78 scenarios recorded from `9077f1e` before any move, compared exactly: results, report texts, exports, search paths, self-tests |
| Science | golden master, primer3-py, ViennaRNA, bcftools, ORFfinder, Mutalyzer and ExPASy comparisons | unchanged reference agreement |
| Components | `tests/test_blast_components.py`, `tests/test_primer_components.py`, `tests/test_primer_screen_configuration.py` | rejected-input messages, executable discovery, local-execution failures and rescue, condition validation, specificity screen, workflow messages, the explicit screen configuration |
| Interface | rendered-menu, dialog and session tests | byte-identical terminal output |
| Packaging | CI import check outside the checkout | every internal package is installed |

The moves themselves were made by a script that cuts a module at
statement boundaries and generates the imports. A second, independent
check found every original statement verbatim in exactly one new place.
The script was not added to the repository: it served this one migration.
The method and the checks are described here and in `docs/CHANGES.md`.

## 6. Remaining architectural limitations

- `analyses/mutation.py` (4,670 lines), `analyses/cloning.py` (3,972),
  `analyses/protein.py`, `analyses/variant_analysis.py` and
  `analyses/sequence.py` are still single modules with built-in
  self-tests. The pattern of this phase applies to them unchanged.
- `LocalBlastDatabaseManager` (985 lines) and `thermodynamics` (1,581
  lines) remain large. They are cohesive, but long.
- Duplicated helpers across modules (`reverse_complement` seven times,
  `normalize_sequence` four, three FASTA readers) are not consolidated
  (audit F6).
- Long dialog functions in `ui/` (audit F7) are unchanged.
- BLAST configuration is still read from the environment at import time
  (audit F9).
- CRISPR code (`analyses/crispr.py`, `ui/crispr_menu.py`) is excluded from
  all of this by owner decision. It does not depend on the BLAST or primer
  modules.
