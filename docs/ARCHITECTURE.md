# Synapse architecture

This describes the code as it is now (last updated 2026-11-20, Phase 10).
It was written in Phase 7.1 as `docs/PHASE71_ARCHITECTURE.md` and is kept
up to date here.
- Findings that led to the changes: `docs/PHASE71_ARCHITECTURE_AUDIT.md`.
- Change records: `docs/CHANGES.md`, sections "Phase 7.1" to "Phase 10".
- Performance on genome-length input: `docs/BENCHMARKS.md`.

## 1. Layers

```
main.py                 launcher from the source folder (python main.py)
synapse_gen/            the terminal session (app.py) and the `synapse` command
ui/                     dialogs and reports (draw through ui/style.py)
core/                   data model, pipeline and provenance, public API, alphabets,
                        storage locations (paths.py), prefix counts for window
                        scans (windows.py)
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

`analyses/blast.py`, `analyses/primer.py`, `analyses/mutation.py` and
`analyses/cloning.py` remain the only modules that callers import:
- the menus in `ui/`, `core/api.py` and `core/pipeline.py`;
- the tests;
- `python -m analyses.<module>`.

Each public module:

- keeps its original docstring, which describes the scientific scope and
  limits;
- re-exports every name it defined before its split, private helpers
  included (`ui/blast_menu.py` uses two of them). Exceptions, each recorded
  in `docs/CHANGES.md`:
  - the removed global `_V5_ACTIVE_CONFIG`;
  - the renamed `_vN_` primer helpers (Phase 8, package B);
- keeps its `if __name__ == "__main__":` block or blocks, and `__all__`
  (`mutation`).

The implementation lives in `analyses/_blast/`, `analyses/_primer/`,
`analyses/_mutation/` and `analyses/_cloning/`.
The leading underscore marks them as internal: their module layout may
change without notice, the public module's names may not.

- **Why two directories.** A package `analyses/blast/` cannot coexist with
  `analyses/blast.py`. Replacing the module would have changed every
  import.
- **Patching in tests.** To patch something a component calls (for
  example `urlopen`, `time.sleep`, `subprocess.run`), patch it in the
  module that calls it:
  - `analyses._blast.ncbi`, `analyses._blast.local_engine`;
  - `analyses._mutation.alignment` (`BANDED_ALIGNMENT_MIN_CELLS`);
  - `analyses._cloning.sites` (`_analyze_enzyme`, `_cached`).

  Patching the public module has no effect on the implementation. Worse,
  a test that does so usually still passes: three tests would have done
  exactly that and were changed (Phase 8).
- **Module state.** Live module state is read where it lives. The cloning
  search cache is `analyses._cloning.cache._STRATEGY_SEARCH_CACHE`; the
  attribute of the same name on `analyses.cloning` is only its value at
  import time.

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
- The historical `_v5`-`_v11` prefixes were replaced by descriptive names
  in Phase 8. The public functions `generate_professional_primer_candidates`
  and `rank_professional_primer_pairs` are defined directly. The "V5"-"V11"
  section comments remain as history.
- The legacy engine (`design`) and the professional engine coexist, as
  before. Both are public API.

## 5. Mutation (`analyses/_mutation/`, Phase 8)

| Module | Lines | Responsibility | Depends on |
|---|---|---|---|
| `models` | 308 | result dataclasses and errors | – |
| `sequence` | 318 | alphabets, IUPAC handling, normalization, translation | `models` |
| `alignment` | 504 | global alignment (banded when exact), identity, edit distance, blocks | `models`, `sequence` |
| `detection` | 396 | mutation calls from an alignment, classification, notation, left normalization | `alignment`, `models`, `sequence` |
| `coding` | 524 | codon and protein consequences, the coding-frame check, protein comparison | `models`, `sequence` |
| `hgvs` | 183 | HGVS c. and p. (3' rule, dup, delins, alleles) | `models` |
| `comparison` | 412 | `detect_mutations`, `compare_sequences`, variant sets, statistics | `alignment`, `coding`, `detection`, `models`, `sequence` |
| `context` | 497 | regions, sequence quality, context, spectrum, hotspots | `models`, `sequence` |
| `impact` | 194 | overlap with restriction sites, primers and guide positions | `models`, `sequence` |
| `reporting` | 523 | validation, text and JSON reports, advanced analysis, map | `context`, `hgvs`, `models` |
| `selftest` | 757 | the built-in self-test | the above |

**Data flow.** `comparison.compare_sequences` -> `detect_mutations`:
1. `alignment.global_align`;
2. `detection.detect_mutations_from_alignment`;
3. `coding.annotate_mutations_with_coding_consequences`;
4. statistics -> `MutationAnalysisResult`;
5. `reporting` and `hgvs` describe the result. `ui/mutation_menu.py`
   draws it.

## 6. Cloning (`analyses/_cloning/`, Phase 8)

| Module | Lines | Responsibility | Depends on |
|---|---|---|---|
| `cache` | 57 | the result cache of one strategy search (deep copies, active only inside `_strategy_search_cache()`) | – |
| `sites` | 446 | sequences, enzyme data, overhangs, recognition sites, restriction maps | `cache` |
| `ends` | 559 | cut ends, their compatibility, junctions | `sites` |
| `strategy` | 1,406 | evaluating and scoring one strategy (digests, internal sites, directionality) | `ends`, `sites` |
| `search` | 830 | searching enzyme combinations, diagnostics, the automatic planner | `cache`, `sites`, `strategy` |
| `summary` | 396 | entry points, summaries, ranking and comparison | `search` |
| `selftest` | 274 | the built-in self-test (checks the API of `analyses.cloning`) | `ends`, `search`, `sites` |

`analyses.cloning` also re-exports three names it imports:
`RESTRICTION_ENZYMES`, `bases_match` and `find_recognition_site_matches`.
The golden master and tests read them there.

## 7. The other analyses (Phase 9)

Same pattern as above: a public module with the original docstring and
entry block, and an internal package that holds the code. Plans:
`scripts/refactor/plans/`.

| Public module | Internal package | Modules, in dependency order |
|---|---|---|
| `analyses/protein.py` | `analyses/_protein/` | `constants`, `sequence`, `translation`, `composition`, `charge`, `profiles`, `report` |
| `analyses/variant_analysis.py` | `analyses/_variant/` | `models`, `vcf`, `normalization`, `annotation`, `statistics`, `export`, `workflow`, `selftest` (`__all__` stays in the public module) |
| `analyses/sequence.py` | `analyses/_sequence/` | `bases`, `composition`, `profiles`, `cpg`, `summary`, `selftest` |
| `analyses/orf.py` | `analyses/_orf/` | `codes`, `features`, `finder`, `summary`, `selection` |
| `analyses/restriction.py` | `analyses/_restriction/` | `sites`, `digest`, `compatibility`, `api`, `selftest` (also re-exports `RESTRICTION_ENZYMES`, `get_restriction_enzyme`, `get_type_iis_enzyme`) |
| `analyses/rna.py` | `analyses/_rna/` | `bases`, `composition`, `coding`, `structure`, `summary`, `pipeline` |
| `analyses/literature.py` (Phase 10) | `analyses/_literature/` | `models`, `normalize`, `http`, `providers`, `ranking`, `enrichment`, `search`, `export`, `context`, `selftest` |

**The literature split (Phase 10).** Phase 9 kept `analyses/literature.py`
whole, because its tests replaced helpers on the module itself. Phase 10
first recorded a characterization section that replaces the network at
`urllib.request.urlopen`, below every Synapse function, then split the
module and pointed every stand-in at the module that calls the helper
(`analyses._literature.providers`, `.search`, `.enrichment`). The
providers' throttle state `_LAST_REQUEST_AT` is a shared dictionary in
`analyses._literature.http`. A stand-in on the public module fails
loudly, because the suite blocks the network.

**Not split.**
- `analyses/crispr.py` is excluded by owner decision.
- Smaller single-file modules stay as they are: `rna_structure`,
  `core/*`, `fileio/*`.

## 8. Testing strategy for refactorings

| Layer | Test | Guards |
|---|---|---|
| Names | `tests/test_module_namespaces.py` | every name of `analyses.blast` (125), `analyses.primer`, `analyses.mutation` (97, and `__all__`) and `analyses.cloning` (61) stays reachable; renamed primer helpers are gone; each internal module of the four packages imports alone and complete; public names are the implementation objects |
| Results | `tests/test_characterization.py` + `tests/fixtures/characterization/` | eleven sections, each recorded before its module moved: `blast_primer.json` (78 scenarios, from `9077f1e`), `mutation.json` (79), `cloning.json` (25), `protein.json` (103), `variant.json` (45), `sequence.json` (111), `orf.json` (27), `restriction.json` (46), `rna.json` (101), `literature.json` (21, from `713a1d4`). Compared exactly, except floats within a relative 1e-12, the last-digit difference of `math.log2` between CPython 3.10 and 3.14 |
| Science | golden master, primer3-py, ViennaRNA, bcftools, ORFfinder, Mutalyzer and ExPASy comparisons | unchanged reference agreement |
| Components | `tests/test_blast_components.py`, `tests/test_primer_components.py`, `tests/test_primer_screen_configuration.py` | rejected-input messages, executable discovery, local-execution failures and rescue, condition validation, specificity screen, workflow messages, the explicit screen configuration |
| Interface | rendered-menu, dialog and session tests | byte-identical terminal output |
| Packaging | CI import check outside the checkout | every internal package is installed |

**How the moves were made.**
- A script cut each module at statement boundaries and generated the
  imports.
- A second, independent check found every original statement verbatim in
  exactly one new place.
- Both scripts and the four split plans are in `scripts/refactor/`, for
  the next modules.

## 9. Remaining architectural limitations

- Large but cohesive:
  - `LocalBlastDatabaseManager` (985 lines);
  - `_primer/thermodynamics` (1,581);
  - `_cloning/strategy` (1,406), whose `evaluate_cloning_strategy` is a
    425-line function.
- Duplicated helpers across modules (`reverse_complement` seven times,
  `normalize_sequence` four, three FASTA readers) are not consolidated
  (audit F6).
- Long dialog functions in `ui/` (audit F7) are unchanged.
- BLAST configuration is still read from the environment at import time
  (audit F9).
- CRISPR code (`analyses/crispr.py`, `ui/crispr_menu.py`) is excluded from
  all of this by owner decision. It does not depend on the split modules.
  `analyses._mutation.impact.crispr_mutation_overlap` was moved unchanged.
- Determinism: all eleven characterization sections are bit-identical under
  different hash seeds (`tests/test_determinism.py` keeps the sequence
  analysis that way). Floats may differ in the last digit between Python
  versions (`math.log2`, compensated `sum()` since 3.12); the snapshot
  comparison allows 1e-12.
