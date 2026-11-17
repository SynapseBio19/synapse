# Phase 6 final assessment

Assessment of Project Synapse as it stands on 2026-11-17, written at the
owner's request before any further development. Analysis and
documentation only: no application code, feature, refactoring or external
validation was changed or run for it.

- Repository: `C:\Users\bayra\OneDrive\Desktop\Project Snyapse - Kopie`,
  branch `main`, HEAD `19f7ef8`, working tree clean.
- Reference specification: the Phase 6 Master Prompt
  (`Project_Synapse_Phase_6_Master_Prompt_Skills_Integrated.txt`), with
  the owner's later decisions: the **AI Experiment Planner stays Coming
  Soon**, and **CRISPR- and CFD-related work stays excluded**. The CRISPR
  module was not examined beyond its existing records.

**Method.**

- Git history read from the start of the completion run (`68e59d1`,
  2026-11-11) to HEAD.
- Every completed work package checked for its code and its tests, and
  each test module run.
- Every row of `docs/FEATURE_STATUS.md` (sections 1-15, 271 rows)
  classified by script and then corrected by hand where the code or the
  tests contradict the row.
- Each Master Prompt section compared with the evidence.
- The full suite run once (516 tests, 0 failures, 5 skipped).
- Earlier measurements reused where none was repeated; their dates are
  given.

---

## 1. Git state and history

| Item | Value |
|---|---|
| HEAD | `19f7ef8` (2026-11-17), working tree clean, no untracked files |
| Commits | 154 in total; 143 in Phase 6 (from `824ca5f`, inclusive) |
| Completion run | `68e59d1`-`ea9612b`, 2026-11-11 to 11-13, 32 commits |
| Continuation run | `224f7d2`-`19f7ef8`, 2026-11-13 to 11-17, 12 commits: 45 files changed, about 3,900 lines added |
| Authors | 150 commits by the owner's Git identity, 4 by "Synapse audit" (the previous cycle); 138 Phase 6 commits carry a `Co-Authored-By` line |
| Integrity | no reset, rebase, force-push or remote push; the original folder `Project Synapse` was only read |
| ReadOnly | 15 of 15 protected files set |

Known record errors (not corrected, because history is not rewritten):

- The test totals in the messages of `70ff2b0` (484) and `a1e3520` (488)
  are one too high; the real totals were 483 and 487. This is recorded
  in `docs/CHANGES.md`.

---

## 2. Are the completed work packages really there?

Each package was checked for its code (a symbol or data entry) and for a
test module that was run now.

| Package | Commit | Code | Tests (run 2026-11-17) |
|---|---|---|---|
| AI Experiment Planner "Coming Soon" | `a2f671b` | `ui/menu.py` entry, `main.py` message | `PlannerPlaceholderTests` 2/2 OK |
| Optional MFE structure (ViennaRNA) | `c5ee3b2` | `analyses/rna_structure.py` | `test_rna_secondary_structure` 12/12 OK |
| bcftools norm comparison; reference-only records reported | `70ff2b0` | `parse_vcf(skipped=)` | `test_bcftools_norm_reference` 7/7 OK |
| ORFfinder comparison | `a1e3520` | fixture with 69 ORFs | `test_orffinder_reference` 4/4 OK |
| HGVS DNA descriptions (Mutalyzer) | `f6f12ed` | `hgvs_c_allele`, `mutation_hgvs_c` | `test_hgvs_mutalyzer_reference` 5/5 OK |
| Instability weights vs ExPASy | `b928fef` | probe fixture | `test_instability_expasy_reference` 3/3 OK |
| `core_nt` | `4a4bb05` | `DATABASES`, catalog | `test_blast_core_nt` 4/4 OK |
| Dialog paths | `c24159f` | – | `test_dialog_paths` 17/17 OK |
| Type IIS enzymes (P2-7) | `4939688` | `data/restriction_enzymes.py` | `test_type_iis` 10/10 OK |
| Circular and 5′-partial ORFs (P1-9) | `f020a1d` | `_origin_spanning_orfs` | `test_orf_partial_circular` 15/15 OK |
| HGVS protein notation | `b7eb1de` | `mutation_hgvs_p` | `test_hgvs_notation` 4/4 OK |
| G·U count (P2-5) | `fedb750` | `hairpin_count_with_gu` | `test_rna_wobble` 3/3 OK |
| Restriction-site primers (10.2) | `99150da` | `design_primers_for_restriction_site` | `test_integration_workflows` 10/10 OK |
| Variant cross-check (P2-4, T12) | `2c68e61`, `cdbf69c` | `cross_check_variants` | `test_variant_crosscheck` 14/14 OK |
| Dialogs moved to `ui/` (P2-2) | `50dea66` | 11 `ui/*_menu.py` | `test_terminal_ui` 34/34 OK (includes `LayeringTests`) |
| BLAST polling resilience (D25) | `a5e5a94` | `await_ready` retries | `test_blast_polling` 5/5 OK |
| Cloning planner cache (P2-6) | `b4b16bb` | search cache | `test_planner_cache` 3/3 OK |
| Primer design speed (P2-6) | `aba3fb0` | screen equivalence | `test_primer_screen_equivalence` 10/10 OK (benchmark opt-in) |
| P0 regressions | Stage 1 | – | `test_p0_regressions` 22/22 OK |
| Golden master | Stage 0 | 12 snapshots | `test_golden_master` 12/12 OK |

**Result:** every package claimed in the final report is present in the
code and covered by passing tests. Exceptions to "passing":

- The four primer3-py comparisons (Tm, hairpin, self-dimer, heterodimer)
  are skipped today, because Windows application control blocks
  primer3-py's compiled module. Their last passing runs are 2026-11-02
  and 2026-11-11.
- The live external checks (ORFfinder, Mutalyzer, ExPASy, NCBI BLAST,
  literature providers) are stored results of single runs, not repeated
  here.

---

## 3. Feature overview

**Unit.** One row of `docs/FEATURE_STATUS.md`: a menu entry, a feature
or a function group, 271 in sections 1-15. The code has 1,151 functions
in `analyses/`, `core/`, `ui/` and `main.py`.

**Categories.**

- **Implemented and tested**: executed by the suite with checked output.
  It also needs one of three things: an independent reference comparison,
  a live check, or the statement that no scientific validation applies.
  Heuristics fall here when they are tested and labeled as heuristics.
- **Implemented but partially validated**: implemented, but only
  rendered or reachable in tests, not tested at all, or without a
  reference where one would apply.
- **Partially implemented**: a documented functional gap.
- **Not implemented**, **Coming Soon** and **Intentionally deferred**:
  as named.

### 3.1 Counts per module

These are the counts after correcting the rows that are out of date
(section 3.3).

| Module | Implemented and tested | Implemented, partially validated | Partially implemented | Not implemented | Coming Soon | Intentionally deferred |
|---|---:|---:|---:|---:|---:|---:|
| Main menu and session | 17 | 0 | 1 (CRISPR entry) | 0 | 1 | 0 |
| Sequence overview | 14 | 0 | 0 | 0 | 0 | 0 |
| ORF analysis | 22 | 0 | 0 | 0 | 0 | 0 |
| RNA analysis | 14 | 0 | 0 | 0 | 0 | 0 |
| Protein analysis | 21 | 0 | 0 | 0 | 0 | 0 |
| Restriction analysis | 18 | 2 | 0 | 0 | 0 | 0 |
| Cloning | 15 | 1 | 0 | 0 | 0 | 0 |
| Primer design | 17 | 2 | 0 | 0 | 0 | 0 |
| Mutation analysis | 20 | 0 | 1 | 0 | 0 | 0 |
| Variant analysis | 21 | 0 | 0 | 0 | 0 | 0 |
| CRISPR (existing records only) | 15 | 0 | 6 | 0 | 0 | 1 (CFD) |
| BLAST | 21 | 3 | 1 | 0 | 0 | 0 |
| Literature | 15 | 1 | 0 | 0 | 0 | 0 |
| Core | 13 | 0 | 0 | 0 | 0 | 0 |
| Support modules | 7 | 0 | 1 | 0 | 0 | 0 |
| **Total (271)** | **250** | **9** | **10** | **0** | **1** | **1** |

### 3.2 Rows not in the first category

**Coming Soon**

- AI Experiment Planner (main menu `[11]`). It states that the feature
  is planned and not yet available, and it runs and stores nothing.

**Intentionally deferred**

- CFD off-target scoring, by owner decision (2026-11-17).
- Genome-wide off-target search belongs here as well. It is recorded as
  a limitation of the CRISPR off-target row and is not a separate row.

**Partially implemented** (10 rows)

| Feature | What is missing |
|---|---|
| CRISPR guide design (menu entry) and six CRISPR rows: ranking, base editing, prime editing, off-target analysis, specificity profile, U6 rule | As recorded before the exclusion: heuristic ranking and editing models; off-target analysis limited to supplied sequences; no CFD. Not re-examined |
| Mutation size guard | Similar long sequences align in an exact band (10 kb in about 0.2 s); unrelated long sequences still need O(n·m) time (about 23 s at 10 kb). There is a note, no limit |
| BLAST provenance | Recorded as partial, but its note lists everything as recorded; whether a gap remains was not settled in this assessment |
| GenBank loading | Feature locations other than spans, `complement()` and `join()` are kept as raw strings, with a note |

**Implemented but partially validated** (9 rows)

| Feature | Why |
|---|---|
| Restriction: complete report; topology selection | Only rendered by the menu tests |
| Cloning: design reasons and orientation view | Only rendered |
| Primer: condition entry; reports | Rendered; the condition prompts are driven with empty answers only |
| BLAST: search configuration | Rendered |
| BLAST: NCBI database download | Not run: a very large download, deliberately avoided |
| BLAST: legacy result file | A labeled pre-correction record; whether to delete it is the owner's decision |
| Literature: Crossref provider | Tested with fixtures only; not part of the live check |

**Not implemented**

- Nothing that is advertised.

**Not implemented and not advertised** (stated as limitations):

- Golden Gate assembly planning.
- Joint normalization of unsplit multi-allelic records.
- HGVS numbering against RefSeq transcripts.
- An organism-mode switch: by design, covered by the genetic code and
  the start policy.
- Genome-scale primer specificity.

### 3.3 Out-of-date rows found in `docs/FEATURE_STATUS.md`

These rows understate the implementation, and the counts above use the
verified state. The document itself was not changed, because this
assessment only documents; updating it is listed in section 6.

| Row | What the document says | Actual state (evidence) |
|---|---|---|
| RNA: Translation & Codon Analysis | "Always NCBI table 1" | Uses the ORF's genetic code since P1-9 (`ui/rna_menu.py` passes `genetic_code`; `test_phase6_tier2`) |
| Cloning: Digest fragments | S17, fragments off by the cut offset | S17 corrected ("The cloning digest cuts where the enzyme cuts", `docs/CHANGES.md`; test in `test_phase6_tier2`) |
| BLAST: Run BLAST Search (local) | Partially implemented | D14-F2 corrected; a real local blastp run is a test (`test_phase6_sequence_types`) |
| BLAST: create custom database | Not tested | `create_custom_database` runs with real `makeblastdb` in `test_phase6_sequence_types` |
| Core: non-interactive API | "Missing: cloning, variant filtering, mutation reporting, BLAST, literature" | All present since S14, plus later entry points (21 in `__all__`) |
| ORF Candidates; mutation coding UI and report; variant manual entry; literature search, configuration and sources | "Reachable" | Driven to the end by `test_dialog_paths` (2026-11-17) |
| Main menu: Literature Search | "Execution-tested (offline)" | Also live-checked on 2026-11-12 |
| Section headers | Line and function counts (for example `primer.py` 5,412 lines) | Out of date (`primer.py` now 5,183 lines) |

---

## 4. Comparison with the Phase 6 Master Prompt

Statuses: COMPLETE, PARTIALLY COMPLETE, NOT COMPLETED, BLOCKED, DEFERRED.
The previous final report's matrix was re-checked.

**One correction:** §12 drops from COMPLETE to **PARTIALLY COMPLETE**,
because the feature matrix has out-of-date rows (section 3.3).

| § | Requirement | Status | Evidence / gap |
|---|---|---|---|
| 0 | Coherent, reliable application; no frontend | COMPLETE | layering tests; no GUI work |
| 1.1, 1.4 | Copy only; preserve work and history | COMPLETE | section 1 |
| 1.2 | BLAST database policy | COMPLETE | `data/blast_db/` excluded; remote live checks; a large local database not obtained (BLOCKED, approval) |
| 1.3, 1.3.1 | External resources traceable; skills and MCP | COMPLETE | ledger with source, date, licence, checksum; 11 skills present in `.claude/skills/`; their use is recorded as consulted, not executed |
| 2 | Full audit, audit report | COMPLETE | `PHASE6_INITIAL_AUDIT.md` |
| 3.1 | P0 re-inspection | COMPLETE | 22 tests |
| 3.2 | ORF | COMPLETE for linear sequences; origin-spanning ORFs not externally compared | ORFfinder 69/69 |
| 3.3 | RNA input rule | COMPLETE | workflow test 10.1 |
| 3.4 | Protein | PARTIALLY COMPLETE | instability weight KQ undecided (ExPASy 24.68, Biopython and Synapse 24.64); primary paper not accessible |
| 3.5 | Primer, Primer3 comparison | COMPLETE | comparison skipped today by application control; passed 2026-11-02/11 |
| 3.6 | Restriction, Type IIS, REBASE | COMPLETE | NEBcutter not run (web application); REBASE and Biopython cover the substance |
| 3.7 | Mutation, HGVS | COMPLETE | Mutalyzer 11/11; numbering is relative to the supplied sequence |
| 3.8 | Variants, bcftools | COMPLETE against bcftools' published test outputs | bcftools itself not run (no Windows build) |
| 3.9 | CRISPR incl. CFD | PARTIALLY COMPLETE; CFD and genome-wide search **DEFERRED** | earlier corrections D8-D10, S3-S5, S23 |
| 3.10 | BLAST | COMPLETE | `core_nt`; live `nt` and `core_nt` |
| 3.11 | Literature | COMPLETE | Crossref not live-checked |
| 3.12-3.13 | Remaining modules; completion procedure | COMPLETE | METHODS per module |
| 4.1-4.3 | Inventory; execution; no disguised placeholders | PARTIALLY COMPLETE | inventory complete but partly out of date (3.3); 9 rows rendered only |
| 5 | Unified terminal UI | COMPLETE (minor remainders) | rendered-menu and dialog tests; self-test "RESULT:" lines in script mode |
| 6 | Output standardization | COMPLETE | Stage G; exports carry parameters and warnings |
| 7.1-7.4 | Core architecture; separation; standards | COMPLETE | `LayeringTests`; pyright 58 (from 134), no new diagnostic; Ruff pyflakes clean |
| 7.5 | No giant modules or functions | PARTIALLY COMPLETE | 4 modules over 3,900 lines; 17 functions of 250 lines or more (largest: `run_mutation_self_test` 689, `main` 632, `find_orfs` 593) |
| 8-9 | Language; style guide | COMPLETE | spelling and prefix rules in the suite |
| 10 | Cross-module workflows | COMPLETE | integration tests 10.1-10.4, 10.7; live 10.5, 10.6 |
| 11.1-11.7, 11.9 | Test levels; golden protection | COMPLETE | 516 tests; golden 12/12 |
| 11.8 | Reproducibility across Python versions | PARTIALLY COMPLETE | 3.10.11 and 3.14.6 run; 3.11-3.13, Linux and macOS not |
| 12 | Feature matrix current | **PARTIALLY COMPLETE** (was COMPLETE) | out-of-date rows, section 3.3 |
| 13-14 | Staged work; boundaries | COMPLETE | one commit per package; no fabrication |
| 15 | Documentation | COMPLETE apart from FEATURE_STATUS | README, METHODS, CHANGES, AUDIT_STATUS, STYLE_GUIDE, TERMINAL_UI_STANDARD exist and match the code where checked |
| 16 | Final criteria | PARTIALLY COMPLETE | section 6 |
| 17 | Final report | COMPLETE | `PHASE6_FINAL_REPORT.md` (A-L) |
| 18 | Final directive | PARTIALLY COMPLETE | open items stated |

---

## 5. Evaluation

### Scientific reliability: good for screening, research and teaching; not validated for decisions

**Agreement with established implementations**, where compared:

- Genetic codes and translation: Biopython.
- Restriction enzyme definitions and digests: REBASE and Biopython.
- Primer Tm: primer3-py, max |ΔTm| 0.0032 °C over 2,916 conditions.
- Protein descriptors: Biopython ProtParam and ExPASy.
- Instability weights: ExPASy, 399/400.
- Variant normalization: bcftools' published outputs, 29/29.
- ORF calls: ORFfinder, 69/69.
- HGVS descriptions: Mutalyzer, 11/11.
- MFE structures: ViennaRNA.

**Weaknesses.**

- Several outputs are heuristic by nature, and are labeled as such:
  primer structure severity, the hairpin screen, guide ranking,
  signal-peptide and transmembrane screens.
- The hairpin screen is sensitive but not specific: only 35 % of its
  candidates are favorable on their own.
- Every comparison used small, chosen datasets: one plasmid, synthetic
  sequences, published test vectors. None measured sensitivity or
  specificity on a large, curated benchmark, such as RefSeq CDS sets for
  ORFs, or GIAB or ClinVar for variants.
- No result is experimentally validated.

### Architecture: clean layering; coarse inside the modules

- Science (`analyses/`, `core/`) runs without the interface, which is
  enforced by tests. The interface goes through one renderer, and there
  is one API for every analysis.
- Inside modules, responsibilities are mixed, and very long functions
  remain (section 4, §7.5).
- There is no separation into packages per domain.
- Optional external tools are integrated cleanly: Biopython, primer3-py
  and ViennaRNA are not runtime dependencies.

### Performance: adequate for single genes and plasmids

| Measurement | Value |
|---|---|
| Primer design on 10 kb | 4.39 s CPU (2026-11-17) |
| Similar 10 kb alignments | about 0.2 s |
| Cloning planner | 3.2 s |
| MFE folding | cubic: 4 kb about 14 s; capped at 3 kb by default |
| Unrelated 10 kb alignments | about 23 s |

Genome-scale work (whole-genome searches, large VCFs) has not been
measured and is not designed for. The pure-Python implementations
would need indexing or native tools for it.

### Test coverage: broad, with known holes

- 516 tests in 36 files: unit, regression, golden master, integration,
  rendered interface, dialog paths, layering, style and external
  reference tests.
- 75 % combined line/branch coverage.
- Weakest: `ui/report.py` 12 % (driven only through a `main.py`
  subprocess, which coverage does not follow), the protein dialog 47 %,
  the BLAST dialog 51 %, `main.py` 34 %.
- The reference comparisons skip when their packages cannot load, as
  primer3-py did today. A run can therefore pass with less evidence than
  usual; the skip reason is printed.

### Maintainability: moderate

**Strengths.**

- Documented conventions, enforced by tests.
- Golden master.
- A change log with reasons.
- Provenance for external data.

**Weaknesses.**

- 52,000 lines with four modules above 3,900 lines.
- 17 functions of 250 lines or more.
- Legacy names (`_v5`-`_v11`, `complexity`).
- 58 remaining pyright errors.
- Documentation spread over about 20 files, which drift (section 3.3).
- No continuous integration.
- No licence.
- A single Windows development machine.

---

## 6. Open requirements and what each needs

| # | Open item | What would close it | Who or what is needed |
|---|---|---|---|
| 1 | FEATURE_STATUS out of date (13 rows, section headers) | Update the rows in 3.3 and the header counts; add a test that checks the header counts against the code | Documentation work, no approval |
| 2 | Instability weight KQ (24.68 vs 24.64) | Read Guruprasad et al. 1990, Protein Eng 4:155-161, Table I; set the published value with a documented golden-master check | Access to the paper (library or the owner) |
| 3 | Origin-spanning ORFs against an external tool | Compare circular calls with ORFfinder's standalone tool or a curated circular plasmid set | Linux or WSL (approval), or a curated dataset |
| 4 | Running bcftools itself | Run the normalization corpus through `bcftools norm` | WSL or Linux (system install, approval) |
| 5 | NEBcutter | Manual digest of pUC19 in NEBcutter V3, compared with Synapse | The owner, by hand; substance already covered |
| 6 | Local BLAST on a real database | Download swissprot (about 0.5 GB) and run blastp | Approval for the download |
| 7 | Genome-scale specificity (primers; CRISPR deferred) | Genome and index, an external aligner, and a design for it | Very large downloads (approval), design work |
| 8 | CFD scoring and other CRISPR items | Lift the owner's exclusion; a citable, licensed matrix | **Owner decision** (deferred) |
| 9 | AI Experiment Planner | A specification of what it should do | **Owner decision** (Coming Soon) |
| 10 | Licence; legacy BLAST result file | Choose a licence; keep or delete `blast_results/TEST_8_REFERENCE_blast.json` | **Owner decision** |
| 11 | Python 3.11-3.13, Linux, macOS | A CI matrix, or local embeddable interpreters for 3.11-3.13 | Small downloads (allowed); CI needs a hosted repository |
| 12 | Large modules and functions; legacy names; `complexity` field | Split by responsibility with byte-identical golden and rendered output; rename with a deprecation path | Refactoring work, no approval |
| 13 | Pyright 58 errors | Triage the remaining Optional and argument-type findings | Work, no approval |
| 14 | Dialog coverage: protein, BLAST, `ui/report.py`, `main.py` | Drive those dialogs like `test_dialog_paths`; run `main.py` in-process for coverage | Work, no approval |
| 15 | primer3-py blocked by application control | Re-run the four comparisons when the block is lifted | Windows policy (not to be bypassed) |
| 16 | Crossref live check; live RID resume | One neutral query; resume a submitted search by RID | Allowed public requests |
| 17 | HGVS against real transcripts | Map supplied CDS to RefSeq transcripts and validate full descriptions with Mutalyzer | Design work; public reference data |
| 18 | Benchmarking on curated data | ORF sensitivity and specificity on a RefSeq CDS set; variant normalization on a GIAB or ClinVar subset; primer success data | Moderate downloads; study design |
| 19 | Experimental validation | Laboratory confirmation of primers and predictions | Outside software |

---

## 7. How far Synapse is from professional, scientifically robust software

**Where it is.**

- Synapse is a carefully engineered research and teaching workbench.
  Its core calculations agree with the reference tools wherever they were
  compared, its heuristics are labeled, and its results are traceable:
  parameters, input hashes and versions.
- It is tested at several levels, and its limits are documented rather
  than hidden.
- For single genes, plasmids and small VCF files it gives answers that a
  trained user can check and use as a first screen.

**What separates it from professional, scientifically robust software**,
in order of weight:

1. **Validation depth.** Agreement on one plasmid and on synthetic cases
   shows that the algorithms are implemented correctly. It does not show
   performance on real data. Professional tools publish benchmarks on
   curated datasets, and most have been reviewed by the community.
   Synapse has neither. This is the largest gap and the slowest to close.
2. **Scale.** Pure Python without indexes is fine at gene and plasmid
   scale, but genomes, large cohorts and genome-wide specificity are out
   of reach. That is acceptable for its purpose, but it excludes many
   professional uses.
3. **Engineering infrastructure.** There is no CI, no testing on several
   platforms or on Python 3.11-3.13, no licence and no release or
   versioning process. Development has run on one Windows machine. These
   are essential for shared or published software and do not need new
   science.
4. **Maintainability debt.** The very large modules and functions, the
   legacy names, the remaining pyright errors and the documentation
   drift raise the cost and risk of every change.
5. **Scope gaps the owner has set aside.** CRISPR scoring (CFD,
   genome-wide off-targets) and the AI Experiment Planner.

**Estimate.**

- *Professional engineering quality* (CI, platforms, licence,
  refactoring, documentation upkeep) is within reach. It would take
  several focused work packages, with no new scientific dependency.
- *Scientific robustness in the professional sense* needs benchmark
  studies on curated data, review by domain experts and, for laboratory
  use, experimental confirmation. That is a project of its own, beyond
  further code changes.

Synapse should not be described as production-ready, clinically
validated, genome-wide validated or equivalent to the tools it was
compared with. It can be described as a tested, transparent platform
whose methods are checked against established tools, suitable for
research support and teaching.

---

## 8. Limits of this assessment

- No new external validation was run. The live results quoted are those
  stored with their dates.
- The CRISPR module was not examined beyond its existing records.
- The feature counts are based on FEATURE_STATUS rows. Only the rows in
  section 3.3 were re-verified against code and tests one by one; the
  other rows were taken as recorded after the classification script, and
  a spot check found no further contradictions.
- Coverage, pyright and benchmark figures are from 2026-11-17.
