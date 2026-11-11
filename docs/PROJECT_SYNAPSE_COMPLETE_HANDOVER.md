# Project Synapse – Complete Handover for a New Session

Written 2026-11-11 at the end of a long Claude Code session, as the primary
context source for the next AI session. It does not replace the
authoritative documents. It explains how they fit together, what is
true now, and what to do next.

**Phase 6 is NOT complete.** Continue it; do not restart it.

## How to read this document

The document separates four kinds of information:

- **[VERIFIED 2026-11-11]**: checked by a fresh inspection when this
  document was written (Git, test runs, pyright, file attributes, file
  sizes).
- **[MEASURED IN SESSION]**: measured earlier in this session and recorded
  in commit messages or documentation, but not re-run for this handover
  (benchmarks, external comparisons).
- **[DOCUMENTED]**: taken from earlier project documents or commit
  messages and not independently re-checked here. This covers the history
  before Phase 6 in particular.
- **[NOT VERIFIABLE]**: stated somewhere but impossible to reconstruct
  from the repository.

Authority order, from highest:

1. The Phase 6 Master Prompt.
2. The source code and tests.
3. `docs/PHASE6_CONTINUATION.md`.
4. This document.
5. All other reports.

A report never overrides the Master Prompt.

## Contents

- Part 0 – Quick start for the next session
- Part I – Project identity and vision
- Part II – Project history and development
- Part III – Modules, one by one
- Part IV – Core architecture and data flow
- Part V – Scientific principles and corrections
- Part VI – Performance work
- Part VII – Testing, validation and quality assurance
- Part VIII – The Phase 6 Master Prompt and the status of each stage
- Part IX – Stage F: terminal UI standardization
- Part X – Development environment and safeguards
- Part XI – Documentation map
- Part XII – Current exact project status
- Part XIII – Task register, blockers and approvals
- Part XIV – Working with the project owner
- Part XV – Exact instructions for the next session
- Part XVI – Completeness check of this document
- Appendix A – Complete commit history
- Appendix B – Source files and sizes
- Appendix C – Public API (`core/api.py`)
- Appendix D – Test files
- Appendix E – Commands
- Appendix F – Status vocabulary and glossary

---

## Part 0 – Quick start for the next session

1. Read Part XV (instructions), Part XII (state) and Part XIII (tasks).
2. Read the Master Prompt:
   `C:\Users\bayra\Downloads\Project_Synapse_Phase_6_Master_Prompt_Skills_Integrated.txt`
   (2,344 lines). It lies **outside** the repository.
3. Read `docs/PHASE6_CONTINUATION.md`.
4. Verify `git log --oneline -1` and `git status`. The HEAD below was
   correct when this document was committed (Part XII).
5. Next task: primer design at 10 kb under 10 s (P2-6). The largest cost is
   `_alignment_stem` in `analyses/primer.py`. See Part XIII, task T1.

---

## Part I – Project identity and vision

### What Synapse is

[VERIFIED from README and code] Synapse is a terminal bioinformatics
workbench in Python 3 for DNA, RNA and protein sequence analysis. Its
scientific core uses only the standard library (`pyproject.toml`:
`dependencies = []`). A user loads FASTA or GenBank records, selects one
sequence, and analyzes it with twelve modules in one session:

- sequence overview;
- ORF search;
- RNA;
- protein;
- restriction and cloning;
- primer design;
- mutation and variant analysis;
- CRISPR guide design;
- BLAST;
- literature search.

Every result is stored in a pipeline with its parameters and an input
hash, and can be exported as a reproducible JSON or Markdown report. A
non-interactive API (`core/api.py`) offers every analysis without
keyboard input.

### Origin and motivation

The repository keeps the origin in `LEGACY/first_dna_analyzer.py`: a
single script that reads a DNA string with `input()`, accepts only A, C,
G and T, and counts nucleotides [VERIFIED]. The project grew from that
script into the current modular application.

The owner's stated intention, in the Master Prompt (section 0) and
throughout this conversation, is to build an intelligent analysis
platform that supports biologists and biotechnology researchers. It
should be a coherent scientific application in which the modules work
together, not a collection of independent scripts.

### Intended users

Biologists, biotechnology students and researchers who need routine
sequence analyses with transparent methods. The README states the
status plainly: "research and teaching software under active
correction. Its results are computational; none has been validated
experimentally" [VERIFIED].

### What distinguishes it

- **Integration.** A selected ORF flows into RNA, protein, CRISPR and
  BLAST. Results are stored with parameters and hashes, and one report
  covers the session.
- **Transparency.** Every heuristic is named as a heuristic, and
  `docs/METHODS.md` gives the equations, defaults, references and the
  validation level of each module.
- **Reproducibility.** Results carry parameters, warnings, the
  enzyme-database version and the input SHA-256. A golden master
  protects behavior.

### What it is not

It does not replace Primer3, NCBI BLAST, ORFfinder, bcftools, REBASE,
NEBcutter, Cas-OFFinder, RNAfold, ProtParam or HGVS validators. Where it
reimplements an established calculation, it compares itself with those
tools where possible and says where it could not.

The AI Experiment Planner is **not implemented**. Its menu entry says so.

It is not a website or GUI project in Phase 6. The Master Prompt forbids
frontend work: §14.1. The folder `Synapse_Website/` exists but is outside
scope.

### Long-term vision

[DOCUMENTED, Master Prompt §0, §7.2, §14.1] A scientifically reliable,
professionally engineered platform whose core is independent of the
terminal interface, so that a later front end (P3-3, deferred) can use
`core/api.py`.

### Philosophy

The Master Prompt's final directive (§18) states the priorities:

- scientific correctness and functional completeness before cosmetics;
- every result traceable;
- every heuristic labeled as such;
- no fabricated data;
- no claims without evidence;
- targeted changes instead of rewrites.

"The priority is not to make Synapse look finished. The priority is to
make Synapse genuinely work."

---

## Part II – Project history and development

The full commit list is in Appendix A. In total there are 125 commits,
from 2026-10-29 to 2026-11-11 [VERIFIED].

### Before Phase 6 (previous cycle, "Claude Cowork") [DOCUMENTED]

According to the Master Prompt §0 and `docs/CHANGES.md`, the previous
cycle produced:

- an audit;
- fixes for the nine P0 defects;
- shared data models;
- pipeline integration;
- `core/api.py`;
- 109 passing tests;
- golden-master snapshots;
- `docs/METHODS.md`.

Git shows this cycle as one day of commits (2026-10-29):

| Commit | Content |
|---|---|
| `e4395ed` | Baseline, the state before the audit-driven changes |
| `f4a2ac0` | Stage 0: test suite, golden master, P0 reproduction tests |
| `08fdce2` | Python 3.10/3.11 import fix in `protein.py` |
| `59d1997` | Stage 1: the nine P0 fixes (see Part V) |
| `b99e63d` | Stage 2: shared data model, visible result storage, file input, logging |
| `5171926` | Stage 3: P1 scientific improvements, reproducible report |
| `dd15af5` | Stage 4: consolidation, shared tables, core API, methods documentation |
| `88cc4ae` | Stage 5: unresolved protein residues keep their position |
| `0dad3bb` | `docs/AUDIT_STATUS.md` |

The detail of that work is known only from these commits and documents
[NOT VERIFIABLE beyond them].

### Phase 6 (2026-10-30 to now)

**2026-10-30.**
- Pyright configuration recorded (`5242fcb`).
- Local skills ignored in Git (`5d8c73a`).
- **Phase 6 initial audit and feature matrix** (`824ca5f`): 19 confirmed
  defects D1-D19 and suspected items S1-S14. Items S15-S24 came later.

**2026-10-31: Stages B, C, D and E work.**
- D1-D19 corrections, one commit each, among them:
  - `dc780f6` D1 initiator Met;
  - `ab4bd2a` D2 an ORF is translated with its own genetic code. This
    commit is also the **pyright baseline: 134 errors**;
  - `d057ac1` central alphabets and the type classifier (D13, D14);
  - `162ca16` one input contract for menu and API.
- Suspected items S1-S24 assessed.
- The superseded primer generations removed (D16 part 2, `dfe773d`).
- RNA codon view with the ORF's genetic code (`8d28586`).
- Scan cache S18 (`87c3fb0`).
- API entry points S14 (`6a8c00c`).
- BLAST RID resume (`ad843c5`).
- `main()` entry point S13 (`7008cd8`).
- Several handover updates.

**2026-11-01: Stage F begins.**
- One separator width (D17, `5db21ab`).
- One meaning of `[0]` and one message prefix (`b1fa927`).
- BLAST menu renumbered (`10431d6`).
- XML exposure measured (S24, `f803577`).

**2026-11-02: Stages G, H and I.**
- No more "confidence" for heuristic scores (D18).
- Evidence levels and unit notation.
- Warnings carried into exports.
- Restriction parameters recorded.
- Editorial pass and `docs/STYLE_GUIDE.md` (Stage H).
- README rewrite.
- Primer Tm compared with primer3-py (`2851e11`).
- IG instability weight corrected from a Biopython reference test
  (`835744e`).
- Workflow tests (`999f77a`).
- AUDIT_STATUS and FEATURE_STATUS updates.
- Python 3.10 f-string fix in CRISPR (`db74504`).
- Final report (`1f0c0f5`) and independent review (`fc77ff4`, `fc13aad`).

**2026-11-02 and 2026-11-03: Stage F completion and its review.** Part IX
gives the detail.
- Re-examination; two conversion rounds covering 16 interface files
  (`2e47ff9` to `f9ed98e`).
- Enter-key regression found by the owner and fixed (`1516339`).
- Closing review and recommendations 1-5 (`46710cd` to `12949f3`):
  - direct prints removed;
  - core no longer imports `ui`;
  - every menu rendered and driven;
  - the single-enzyme cloning crash fixed (`8eddc72`; it had existed
    since the audit commit).

**2026-11-11: autonomous completion run, packages A and B** (`68e59d1` to
`690c9e8`). Part XII describes it; the owner stopped it after package B.

---

## Part III – Modules, one by one

The line counts are [VERIFIED 2026-11-11]. "Status" refers to
`docs/FEATURE_STATUS.md` unless stated otherwise. Note that FEATURE_STATUS
has **not yet** been updated for the commits of 2026-11-11 (Part XI).

### Common model

- **Input.** `fileio/loader.load_sequences()` reads FASTA (`fileio/fasta.py`)
  or GenBank (`fileio/genbank.py`). `core/alphabets.classify_sequence()`
  assigns DNA, RNA or PROTEIN with an evidence level: CONFIRMED,
  AMBIGUOUS or INVALID. An ambiguous record is declared by the user
  (`main.resolve_sequence_type`).
- **Records.** `core/models.SequenceRecord(identifier, sequence, description,
  topology, alphabet, ...)`.
- **Results.** `core/models.AnalysisResult` (module, data, parameters,
  warnings, limitations, input SHA-256). It is produced by `core/api._result()`
  and stored by `core/pipeline.PipelineContext.store()`.
- **Interface.** Every screen uses `ui/style.py`; Part IX describes the
  renderers.

### 1. Sequence overview – `analyses/sequence.py` (2,493 lines), `ui/report.py` (721)

- **Purpose.** Composition, GC/AT content and skew, ambiguity codes,
  homopolymers, CpG islands, profiles, molecular mass and "complexity".
- **API.** `core.api.sequence_overview()`. The screen is
  `ui.report.print_results`.
- **Science.** "Complexity" is the normalized Shannon entropy of the A/C/G/T
  composition. It does **not** detect repeats: `ACGT`×100 scores 1.000.
  METHODS states this (P2-5).
- **Status.** Implemented and execution-tested.
- **Known limitation.** The field name "complexity" was kept for
  compatibility.

### 2. ORF analysis – `analyses/orf.py` (3,360)

- **Purpose.** ORF search on both strands; candidate comparison,
  relationships (overlap and nesting), selection, settings.
- **Key functions.**
  - `find_orfs(dna, include_partial=True, min_length=75, start_codons, stop_codons, start_codon_policy, nested=False, ranking="length", genetic_code=1, topology="linear", include_5prime_partial=False)`;
  - `analyze_orfs`;
  - `analyze_orfs_enhanced`;
  - `select_orf`;
  - `translate_orf_sequence`;
  - `orf_has_start_codon`;
  - `orf_end_on_sequence`;
  - `_keep_longest_per_stop`;
  - `_five_prime_partial_orfs`;
  - `_origin_spanning_orfs`;
  - the menu `run_orf_analysis_ui(..., topology)` and its settings
    `_orf_settings_menu`.
- **Defaults** (ORFfinder conventions): 75 nt including the stop, ATG only,
  the longest ORF per stop, ranking by length. `LEGACY_ORF_SETTINGS`
  reproduces the behavior before 2026-10-29.
- **Genetic codes.** NCBI tables 1, 2, 4 and 11 via `core/seqtools.py`. An ORF
  dict carries its `genetic_code`, and downstream translation uses it
  (D2).
- **Coordinates.** 0-based inclusive `start <= stop` on the forward strand
  for both strands.
- **New on 2026-11-11 (P1-9).**
  - `topology="circular"` adds ORFs across the origin (`wraps_origin`,
    unrolled `stop`, `stop_on_circle`) and drops partial ORFs.
  - `include_5prime_partial` adds "Partial - Missing Start" and
    "Partial - Missing Start and Stop" records. Their first codon is never
    read as initiator Met.
  - The per-stop key no longer contains the frame.
- **Data exchange.** The selected ORF dict (`sequence`, `strand`, `frame`,
  `start`, `stop`, `genetic_code`, `status`, `search_parameters`) goes to
  RNA, protein, CRISPR and BLAST. `core/orf_context.resolve_orf_region()`
  validates its coordinates against the parent sequence.
- **Tests.**
  - golden master (`orf.json`, `orf_legacy_baseline.json`);
  - `test_orf_partial_circular.py`, including rotation invariance;
  - tier 1 and 2 tests;
  - comparison with Biopython for the genetic codes.
- **Not done.** No ORFfinder comparison (blocked). Organism mode is
  covered by genetic code plus start policy, with no separate switch.

### 3. RNA analysis – `analyses/rna.py` (3,093)

- **Purpose.** Transcription; an overview with GC/AU content; a
  three-frame codon view; features (hairpins, homopolymers, poly-A,
  start/stop landscape); a report; a screening score.
- **Rule (Master Prompt §3.3).** With a selected ORF the input is that ORF.
  Without one it is the whole selected DNA sequence. The menu shows which
  source it uses ("Input source").
- **Codon view.** It uses the selected ORF's genetic code (P1-9 part,
  `8d28586`).
- **Hairpins.** `find_rna_hairpins(rna, min_stem=4, min_loop=3, max_loop=12, allow_gu=False)`
  screens stems of 4-7 nt that pair at 75 % or more, Watson-Crick only.
  Since `fedb750`, `summarize_rna_structure()` adds
  `hairpin_count_with_gu`, and the screen shows "Including G-U pairs"
  (58 → 260 on the golden reference).
- **API.** `core.api.analyze_rna`.
- **Tests.** Golden master (`hairpin_count` 58) and `test_rna_wobble.py`.
- **Limitation.** No free-energy model and no RNAfold (an external binary,
  not installed).

### 4. Protein analysis – `analyses/protein.py` (4,369)

- **Purpose.** Translation and physicochemical descriptors: mass, pI,
  GRAVY, instability index, aliphatic index, aromaticity, extinction
  coefficient, composition, hydropathy and flexibility windows, and
  transmembrane and signal-peptide screens.
- **Entry points.**
  - `run_protein_analysis_from_orf(orf_sequence, orf_id, strand, reading_frame, orf_status, pipeline, genetic_code)`;
  - `run_protein_analysis_from_record` for protein records;
  - `analyze_complete_protein`;
  - `analyze_translation(rna, initiator_start, genetic_code)`;
  - API: `analyze_protein_from_cds`, `analyze_protein`.
- **Rules.**
  - The initiator Met applies only to an initiation codon of the code (D1).
  - Unresolved residues (`?`, B/Z/J/X) keep their position (Stage 5,
    D14-F6).
  - U and O are defined residues.
- **Validation.** Compared with Biopython ProtParam on 31 proteins and
  ExPASy on ubiquitin. The IG instability weight was corrected from
  −7.49 to 1.0 against Biopython's DIWV table.
- **Limitation.** The Guruprasad 1990 primary table was not consulted
  (blocked).
- **Since 2026-11-11.** A 5′-partial ORF is translated without an
  initiator.

### 5. Restriction analysis – `analyses/restriction.py` (2,830), `data/restriction_enzymes.py` (421)

- **Enzyme database.** 22 palindromic enzymes in `RESTRICTION_ENZYMES`
  (`ENZYME_DATABASE_VERSION = "synapse-2026-10-29"`).
- **Type IIS (`4939688`).** `TYPE_IIS_ENZYMES` (BsaI, BsmBI, BbsI, SapI;
  `TYPE_IIS_DATABASE_VERSION = "synapse-2026-11-11"`) is a separate table.
- **Matcher.** `find_recognition_site_matches()` with the cached
  `_scan_recognition_sites()`. A site is CONFIRMED (exact bases) or
  POSSIBLE (IUPAC ambiguity); only CONFIRMED sites cut (P0-5).
- **Main functions.**
  - `perform_digest(sequence, enzyme_names, topology, include_possible)`;
  - `get_cut_positions`;
  - `find_type_iis_sites`;
  - `create_restriction_map`;
  - `analyze_restriction_site`;
  - `calculate_gel_bands`;
  - unique and non-cutters;
  - `recommend_enzyme_pairs`;
  - `analyze_restriction`.
- **Menus.** `run_restriction_cloning_ui` → `run_restriction_analysis_ui`
  (overview, map, digest and gel, complete report).
- **Topology.** Linear and circular; origin-spanning sites (D3).
- **Coordinates (`1ae2833`).** In site and cut records, unsuffixed
  `site_position`, `top_cut` and `bottom_cut` are 1-based (for a cut, the
  first base after it). `*_0based` keys give the 0-based start or the cut
  boundary. `cut_positions` lists boundaries.
- **Validation.** Biopython 1.88 `Bio.Restriction` (REBASE-derived,
  local data):
  - definitions 22/22;
  - pUC19 digests 132/132;
  - end compatibility 484/484;
  - Type IIS definitions 4/4;
  - 32 digests with 140 cuts and 0 differences.

  It was not compared with NEBcutter (blocked).
- **Limitation.** Methylation is shown as a note for ClaI, XbaI, DpnI and
  HpaII; blocking and star activity are not modeled.

### 6. Cloning – `analyses/cloning.py` (5,057)

- **Purpose.** Single-enzyme analysis, two-enzyme strategies and an
  automatic planner over all enzyme combinations.
- **Functions.**
  - `analyze_cloning`;
  - `evaluate_cloning_strategy`;
  - `automatic_cloning_planner`;
  - `find_best_cloning_strategies`;
  - `ends_compatible`;
  - `print_cloning_report` with `_print_single_enzyme_report` (`8eddc72`)
    and `_print_single_strategy_report`.
- **Science.** End compatibility requires the same overhang sequence and
  the same polarity (5′ vs 3′); blunt ends are compatible with blunt ends
  (P0-3). Ambiguous sites never define ends (D4). Vector and insert
  topology are respected (D5).
- **Performance (`b4b16bb`).** A search-scoped cache
  `_strategy_search_cache()` with `_cached()` and `_copy_result()`. The
  planner went from 21.7 s to 3.2 s with byte-identical results.
- **Type IIS.** Deliberately excluded from the fixed-overhang model.
  Golden Gate planning is not implemented.
- **Known text issue.** The score reasons say "two-end architectures" even
  in single-enzyme mode. This was left unchanged and is noted.

### 7. Primer design – `analyses/primer.py` (5,058)

- **Purpose.** Primer pairs for a target region, single-primer analysis
  and the Tm-model dialog.
- **API.** `core.api.design_primers` → `design_primers_professional`, which
  returns a `ProfessionalPrimerDesignReport` with `ranked_pairs` of
  `AdvancedPrimerPair(pair=PrimerPair(...))`.
- **Thermodynamics.**
  - SantaLucia 1998 unified nearest-neighbor parameters;
  - salt correction 0.368·(N−1)·ln[Na⁺eq];
  - Mg²⁺ after von Ahsen 2001: Na⁺eq = Na⁺ + 120·√(Mg²⁺ − dNTP);
  - defaults 50 mM Na⁺, 1.5 mM Mg²⁺, 0.2 mM dNTP, 250 nM;
  - alternative parameter sets are kept.
- **Validation.**
  - The Primer3 manual example is reproduced.
  - Compared with primer3-py 2.4.0 on 2,916 conditions: maximum |ΔTm|
    0.0032 °C.
- **Structures.** `analyze_hairpin`, `analyze_self_dimer` and
  `analyze_heterodimer` classify the longest perfect stem: ≥ 4 / 7 / 10 bp
  is MODERATE / HIGH / VERY HIGH, together with 3′ complementarity
  (`_v8_structure_severity`). This is a **heuristic**. Against primer3-py
  ΔG the Spearman ρ is −0.78 (hairpin), −0.74 (self-dimer) and −0.51
  (heterodimer), and the classes overlap (`27255df`,
  `tests/test_primer_structure_reference.py`).
- **Coordinates.** `PrimerCandidate.start` and `end` are 0-based half-open
  on the template. The menu shows 1-based positions.
- **Since 2026-11-11.**
  - `primer_pair_binding_sites()` (P2-4).
  - The hairpin search slices one reverse complement and stops at the
    first match (`15e2ff1`): 10 kb design went from 135.6 s to 25.6 s.
- **Open.** The 10 kb design is still above the 10 s target (T1). The
  v5–v11 ranking chain is kept (P2-1 partial). There is no
  restriction-site-driven primer design (T2).

### 8. Mutation analysis – `analyses/mutation.py` (5,273)

- **Purpose.** Compare a reference and a variant sequence: global
  alignment, mutation list, statistics, coding consequences, report.
- **Functions.**
  - `compare_sequences`;
  - `detect_mutations`;
  - `global_align`;
  - `_banded_align`;
  - `_banded_is_exact`;
  - `sequence_similarity`;
  - `_edit_distance`;
  - `analyze_coding_consequences`;
  - `mutation_hgvs_c` and `mutation_hgvs_p`;
  - `format_mutation`;
  - `build_mutation_report`;
  - the menus `run_mutation_variant_ui` and `run_mutation_analysis_ui`.
- **Alignment.** Needleman–Wunsch with match 2, mismatch −1 and gap −2.
  Ties are broken diagonal > up > left. Since `f9d4356` it keeps two score
  rows and a one-byte traceback per cell. Since `e5bb456`, above 1,000,000
  cells it uses an exact banded alignment with an optimality bound
  (Part VI).
- **Similarity.** Levenshtein distance, exactly banded above 1,000,000
  cells (`15e2ff1`).
- **Coding consequences.** Position 1 is assumed to be codon 1, and the
  CDS can come from the selected ORF. Start-loss is reported only for
  codon 1 (P0-6).
- **HGVS (`b7eb1de`).** Predicted changes are written `p.(R12C)`,
  `p.(R12=)`, `p.(W24*)`, `p.(M1?)` and `p.(*110Qext*?)`; they appear in
  report lines. This is HGVS-style, not validated with a validator.
- **Navigation (`8379ac2`).** `_prompt_sequence` offers `[0] Cancel`, which
  raises `SequencePromptCancelled`.

### 9. Variant analysis – `analyses/variant_analysis.py` (2,910)

- **Purpose.**
  - VCF 4.x parsing (text, gzip, BGZF);
  - INFO and FORMAT fields, cardinality Number=A/R/G, multi-allelic
    splitting, genotypes and phasing;
  - FILTER semantics (`.` = not applied);
  - REF validation (`check_reference_alleles`);
  - left-normalization (`left_normalize_variant`, after Tan et al. 2015);
  - filtering, statistics (Ti/Tv, allele frequencies, quality), cohort
    summary, export (VCF, TSV, report).
- **Menu.** It follows the Master Prompt §5.3 order: input 1–2, normalize 3,
  filter 4, statistics 5, analyze 6, compare 7, cohort 8, export 9–11,
  self-test 12.
- **Contig matching.** Exact by default; `chr_prefix` is opt-in. A missing
  reference is flagged (D12).
- **Bridge.** `annotate_from_mutation_analysis()` stores the class in
  `consequence` and, since `b7eb1de`, the `p.(...)` notation in
  `protein_change`.
- **P2-4 (`2c68e61`).** `check_features_against_variants()` compares
  variants with primer and guide sites. The API entry point is
  `core.api.cross_check_variants`. It is **not** wired into the menus,
  because the pipeline does not record which VCF contig was mapped to the
  loaded sequence.
- **Not done.** No `bcftools norm` comparison (blocked).

### 10. CRISPR guide design – `analyses/crispr.py` (3,083)

- **Systems.** SpCas9 (NGG), SpCas9-NG (NG), SaCas9 (NNGRRT) and Cas12a
  (TTTV, 5′ PAM).
- **Functions.**
  - `design_guides`;
  - `run_crispr_analysis`;
  - `rank_guides`;
  - `pam_proximal_seed`;
  - `cut_site_forward`;
  - `search_local_off_targets`;
  - `specificity_profile`;
  - base editing (`analyze_base_editing`; CBE/ABE windows);
  - prime editing (`design_prime_editing_candidates`);
  - exports (JSON, CSV, FASTA, batch);
  - `guide_binding_site` (P2-4).
- **Menu.** `crispr_menu` with sequence and region, batch, system, target
  discovery, guide analysis, editing, specificity and reports.
- **Corrections.**
  - The modeled cut site is reported with its source (D8).
  - A guide's own site is excluded from off-targets (D9).
  - The configuration fields are honored (D10).
  - The seed is PAM-proximal for Cas12a (S23).
  - The poly-T (U6) flag is explained (S5).
- **Coordinates.** `guide_start`, `guide_end`, `pam_start` and `pam_end` are
  1-based inclusive.
- **Limitations.**
  - The off-target search covers only supplied sequences, not a genome
    (P3-2 deferred).
  - CFD scoring is rejected without the published matrix (P2-3, blocked).
  - The ranking is a transparent heuristic.

### 11. BLAST – `analyses/blast.py` (5,576)

- **Engines.** Local BLAST+ (`LocalBlastEngine`, BLAST+ 2.17.0 installed) or
  remote NCBI (`NCBIBlastClient`: submit/status/fetch/run/resume). The
  remote client asks for confirmation and shows a privacy notice.
- **Database manager.** `blast_database_manager_menu`,
  `LocalBlastDatabaseManager`; the folder `data/blast_db/` is git-ignored
  and deliberately absent.
- **Parsing.** `parse_blast_xml` is tested independently of execution.
- **Query coverage.** The union of HSP query intervals divided by the query
  length (P0-1).
- **Evidence profile.** Bands are gated at E ≤ 1e-3; an insignificant hit
  is not treated as evidence.
- **Since 2026-11-11.** HTTPError responses are closed (`68e59d1`).
- **Not done.** No live remote run in Phase 6 (approval needed). The file
  `blast_results/TEST_8_REFERENCE_blast.json` is a labeled legacy record
  (D7); deleting it is the owner's decision.

### 12. Literature search – `analyses/literature.py` (2,395)

- **Providers.** PubMed, Europe PMC, OpenAlex and Crossref, with
  deduplication.
- **Requests.** Per-host throttling and bounded retry on 429/5xx with
  `Retry-After` (`_request_text`, `_retry_delay`). The HTTPError response
  is closed (`68e59d1`).
- **Data.** `LiteratureQuery` and `LiteratureSearchRun`; the date
  validation flags future dates. The menu is `_interactive_menu`.
- **Rule.** No fabricated publications.
- **Status.** The offline tests pass. Live providers have not been run (S11
  remains unverified until the owner approves).

### 13. AI Experiment Planner

Main menu `[11]` is labeled "not implemented" (D15). P3-1 is deferred. No
code exists for it [VERIFIED].

### 14. Pipeline – `core/pipeline.py` (1,317)

- **Classes.** `PipelineContext` and `SynapsePipeline`.
- **Run methods.** `run_sequence`, `run_orf` (passes the menu settings and,
  since 2026-11-11, `topology`), `run_rna`, `run_protein`,
  `run_restriction`, `run_cloning`, `run_cloning_planner`.
- **Selection.** `set_selected_orf` / `select_orf`.
- **Storage.** `store(module, data, parameters, warnings)` keeps
  `AnalysisResult` objects per module.
- **Status.** `status_report()` returns data, and `print_status(render)`
  hands it to a `StatusRenderer` (`b25e800`). The terminal renderer is
  `ui.report.print_pipeline_status`.
- **Rules.** Selecting an ORF changes the sources of RNA and protein. Stale
  dependent results are invalidated (`ORF_DEPENDENT_MODULES`). The protein
  step uses the selected ORF's genetic code; since 2026-11-11 a 5′-partial
  ORF is translated without an initiator.
- **Report.** `core/report.py` builds and exports JSON or Markdown
  (`build_analysis_report`, `export_analysis_report`, `canonical_json`).

### 15. Other components

| File | Content |
|---|---|
| `core/alphabets.py` | Classifier and IUPAC sets |
| `core/seqtools.py` | Genetic codes, reverse complement, GC |
| `core/orf_context.py` | ORF region resolution |
| `core/models.py` | Data model |
| `core/validation.py` | Legacy DNA validation |
| `core/logging_setup.py` | Logging setup (`SYNAPSE_LOG_DIR`) |
| `data/amino_acids.py` | Amino-acid data |
| `data/genetic_code.py` | Genetic-code data (4 lines) |
| `ui/menu.py` | `show_menu` and `select_gene` |
| `ui/report.py` | Sequence overview screen and pipeline status |
| `ui/style.py` | All shared renderers |

---

## Part IV – Core architecture and data flow

### Layers [VERIFIED]

```
main.py                session loop: load -> type gate -> select -> main menu
ui/style.py            renderers (header, section, field, menu, ask, confirm, messages, columns)
ui/menu.py, ui/report.py
analyses/*.py          science + (still) their own menus  <- P2-2 open
core/api.py            non-interactive entry points, AnalysisResult
core/pipeline.py       PipelineContext / SynapsePipeline, result storage, status as data
core/models.py, core/alphabets.py, core/seqtools.py, core/orf_context.py, core/report.py
fileio/                FASTA / GenBank reading
data/                  enzymes, genetic codes, amino acids
```

### Entry point

`main.py:main()` (line 305; guarded by `if __name__ == "__main__"`, S13). The
session runs as follows:

1. `choose_sequence_file()`: Enter loads `examples/dna.fasta`, Q quits, and
   a missing file is reported.
2. `load_sequences()`.
3. `select_gene()`.
4. `resolve_sequence_type()`.
5. For a protein record, `run_protein_record()`; otherwise the main menu.
   It has entries 1–11 and the keys P, R, S and 0.

### Dependency rules [VERIFIED by tests]

- `ui/style.py` imports nothing from `analyses/`, `core/`, `data/` or
  `fileio/`.
- `core/` imports `ui` nowhere (`LayeringTests`, including a subprocess
  that refuses `ui`).
- **Open (P2-2).** All twelve analysis modules import `ui.style` at module
  level because their menus live in them: 133 top-level definitions, about
  8,900 lines. Measured in this session. `core` itself runs without `ui`,
  but an analysis called through `core.api` imports `ui` with its module.

### Data-flow rules

- **RNA source.** It is the selected ORF if one is selected, otherwise the
  whole DNA. Never switch silently; the source is shown and recorded
  (Master Prompt §3.3).
- **Genetic code.** It travels with the ORF dict (`genetic_code`) into
  `translate_orf_sequence`, `analyze_protein_from_cds`, the pipeline
  protein step and the RNA codon view. An explicit conflicting table
  raises an error (D2).
- **Sequence type.** The classifier's result (`alphabet`,
  `type_evidence`, `declared_type`) travels with `SequenceRecord`. The API
  refuses wrong types through `_sequence_of(..., accepted=...)`.
- **Topology.** It travels with the record into restriction, cloning and,
  since 2026-11-11, ORF search (pipeline, API, menu).
- **ORF coordinates.** Validated against the parent before CRISPR and
  BLAST use them. An origin-spanning ORF falls back to ORF-relative
  positions with a note.

### Error handling

- Modules raise their own error classes (`MutationAnalysisError`,
  `CRISPRAnalysisError`, `BlastAnalysisError`, `VariantAnalysisError`,
  `LiteratureSourceError`, ...).
- Menus print `ERROR:` messages through `style.error` and stay usable.
- `main.report_error()` frames unexpected errors.

### External services

NCBI BLAST, the four literature providers and `update_blastdb.pl`
downloads. None runs without explicit user action in the menus. In
development they follow the approval policy in Part X (owner decision of
2026-11-12): neutral requests are allowed, sending private sequences or
data needs approval.

### Implemented architecture improvements

- one input contract (D13/D14);
- central alphabets;
- one restriction matcher;
- the API for every analysis (S14);
- `main()` entry (S13);
- the shared UI module;
- status as data (core free of `ui`);
- a search cache in cloning;
- one coordinate convention in restriction.

### Outstanding architecture work

- **P2-2.** Move the menus out of the analysis modules.
- **P2-1.** Merge the primer generations.
- Pipeline menu-path parameters for sequence overview and cloning.
- Wiring P2-4 into the menus. This needs the VCF contig mapping stored.

---

## Part V – Scientific principles and corrections

### Principles in use

- **Alphabets.**
  - DNA: A, C, G and T, plus the IUPAC codes R Y S W K M B D H V N.
  - RNA: the same with U.
  - Protein: 20 standard residues, U (Sec), O (Pyl), the ambiguity codes
    B/Z/J/X, `?` (unresolved) and a terminal `*`.
  - A sequence whose letters fit several alphabets is AMBIGUOUS and must
    be declared. It is never guessed (D13/D14).
- **Ambiguity in matching.**
  - Restriction sites: an IUPAC symbol inside a site gives POSSIBLE, never
    a cut.
  - Translation: `?` unless every expansion encodes one amino acid.
  - Alignment: IUPAC-aware match (N matches any base).
- **Genetic codes.** NCBI tables 1, 2, 4 and 11.
  - Initiation codons per table; table 1 adds Synapse's GTG policy.
  - The first codon is read as Met only if it is an initiation codon of
    the table (D1). Internal codons use the elongation code.
- **ORFs.** Start codon to in-frame stop, inclusive. The defaults are
  screening conventions (ORFfinder-like), not biological thresholds.
  - 3′-partial ORFs (no stop) are reported by default.
  - 5′-partial ORFs are optional.
  - On a circle, ORFs across the origin are found by searching
    `dna + dna[:n-1]`.
- **Transcription.** T → U of the coding strand. The input is the selected
  ORF or the whole sequence (Part IV).
- **Restriction.**
  - `cut_top` and `cut_bottom` are offsets from the site start in forward
    orientation. `cut_top < cut_bottom` means a 5′ overhang;
    `cut_top > cut_bottom` means a 3′ overhang.
  - Type IIS enzymes cut outside their asymmetric site. They are searched
    in both orientations: for a forward site at p the cuts are
    p + cut_top and p + cut_bottom; for a reverse site they are mirrored.
    The overhang depends on the sequence.
  - End compatibility requires the same overhang sequence and polarity.
- **Primer thermodynamics.** See Part III.7. The predicted Tm is a model
  value, not an annealing temperature. The structure severities are stem
  heuristics.
- **Alignment.** Global Needleman–Wunsch with a linear gap cost. The
  similarity is 100·(1 − edit distance / max length).
- **HGVS.** HGVS-style notation (den Dunnen et al. 2016): `c.` with
  positions in the supplied sequence; predicted protein changes in
  parentheses. Not checked with a validator.
- **Variants.** VCF positions are 1-based. Left-normalization follows
  Tan et al. 2015 and was checked by hand on constructed cases.
- **CRISPR.**
  - PAM recognition per system.
  - The cut is modeled from `system.cut_after` (blunt for Cas9, staggered
    for Cas12a).
  - The seed is the PAM-proximal part of the guide.
  - Base-editing windows are protospacer positions.
  - The off-target search is local only.
- **RNA hairpins.** A Watson-Crick stem screen without energies. G·U wobble
  pairs are counted separately.
- **BLAST.** E-values depend on query length and database size (METHODS);
  only significant hits get evidence bands.
- **Literature.** Results are retrieved metadata. Relevance is text
  overlap, not study quality.

### Corrections (what was wrong, what changed, how validated, what remains)

**P0 (previous cycle, Stage 1 `59d1997`) [DOCUMENTED].** The nine P0
defects of the Master Prompt §3.1:

1. BLAST query coverage;
2. CRISPR selected-ORF coordinate mapping;
3. cloning end compatibility and overhang polarity;
4. enzyme definitions;
5. ambiguous recognition sites;
6. start-loss interpretation;
7. variant normalization and REF validation;
8. alternative start translation;
9. ORF minimum length and ranking.

`tests/test_p0_regressions.py` keeps 22 known-answer tests for them.

**Phase 6 D items** (`docs/PHASE6_FINAL_REPORT.md` §B, where every item
has its evidence):

| Item | Was wrong | Correction |
|---|---|---|
| D1 | Initiator Met for any GTG/TTG start, whatever the genetic code | Only initiation codons of the code |
| D2 | Translation fell back to table 1 | An ORF carries its own code |
| D3 | Map and sites ignored circular topology | Topology honored |
| D4 | Ambiguous sites could define cloning ends | They cannot |
| D5 | No circular vector or insert | Supported |
| D6/S7 | A local BLAST database given by path always failed | Fixed |
| D7 | A wrong BLAST result file shipped unlabeled | Labeled as legacy |
| D8 | No CRISPR cut site anywhere | Modeled and reported |
| D9 | A guide's own site counted as an off-target | Excluded |
| D10 | Five CRISPR configuration fields were never read | Honored |
| D12 | Contig mismatches were skipped silently | Flagged |
| D13/D14 | A peptide written in IUPAC letters loaded as DNA; two `normalize_dna` contracts | Classifier, one contract |
| D14-F6 | Protein residues B/Z/J/X/U/O lost or shifted | Kept in place |
| D15 | Planner entry looked implemented | Labeled "not implemented" |
| D16 | Shadowed primer definitions | Removed after a reachability review; bytecode-identical survivors |
| D17 | Separators of 11 widths; 7 meanings of `[0]` | One width; one meaning |
| D18 | Heuristic scores called "confidence" | `heuristic_score` (golden key renamed, values proven equal) |
| D19 | Pipeline ORF runner ignored the menu's search settings | Fixed |

The S items (S3–S24) are in the same table. Among them:

- the variant-to-mutation bridge never worked (S15/S16);
- the cloning digest cut at the site start (S17);
- the start classification used a fixed codon set (S19);
- the CRISPR seed end was wrong for Cas12a (S23).

Further corrections:

- **IG instability weight.** It was −7.49; it is now 1.0, the value in
  Biopython's DIWV table. Found by a reference test; the primary paper
  was not consulted.
- **Python 3.10.** An f-string made the CRISPR module unimportable on 3.10
  and 3.11. Found by a static check; there has been no runtime 3.10 test.

**2026-11-11 and the Stage F review (2026-11-03).** Tests for these are in
the suite [VERIFIED]; the measurements are [MEASURED IN SESSION].

| Commit | Was wrong | What changed and how it was validated | Limitation |
|---|---|---|---|
| `8eddc72` | Single-Enzyme Analysis crashed (`KeyError: 'left'`) since the audit | Own report from existing fields | – |
| `1516339`, `c78800e` | Stage F rewriting broke Enter defaults in 15 prompts and the Tm dialog | Restored; `EmptyAnswerTests` | – |
| `68e59d1` | HTTPError responses (literature, BLAST) were never closed | Closed; the tests fail on the old code | – |
| `8379ac2` | The mutation prompt had no `[0]`; ten errors were drawn as fields | `[0] Cancel`; ERROR lines; a rule finds all 10 old sites | – |
| `4939688` | No Type IIS enzymes; forward-only search | Both strands; Biopython: 140 cuts, 0 differences | No Golden Gate planning; not in cloning |
| `b7eb1de` | `p.R2C`, "Ter" mixed into one-letter code; helper unused; variant `protein_change` held the class | `p.(R2C)` etc.; used in reports and variants | Not validated with an HGVS validator |
| `1ae2833` | `site_position` 1-based in maps, 0-based in digests | Unsuffixed keys are 1-based everywhere | – |
| `fedb750` | G·U pairs ignored in RNA hairpins | Separate count (58 → 260); screen unchanged | No energy model |
| `f020a1d` | No origin-spanning or 5′-partial ORFs | Implemented; rotation-invariance test | No ORFfinder comparison |
| `2c68e61` | No primer or guide checks against variants | API-level check | Not in the menus |

### Established, heuristic, predicted

| Kind | Examples |
|---|---|
| **Established calculations** | Translation with NCBI tables; restriction cut positions; nearest-neighbor Tm; Needleman–Wunsch; Levenshtein distance; pI and molecular-weight formulas |
| **Heuristics** | ORF heuristic score; primer structure severity; CRISPR ranking; RNA hairpin screen; transmembrane and signal-peptide screens; BLAST evidence bands; literature relevance; "complexity" |
| **Predictions and approximations** | Tm (not annealing temperature); modeled CRISPR cut; HGVS consequences; protein descriptors |

Nothing is validated experimentally.

---

## Part VI – Performance work

All figures are [MEASURED IN SESSION] on the owner's Windows machine,
recorded in the commit messages and METHODS, and not re-run for this
handover. Results stayed identical in every case.

| Item | Before | After | Evidence of identical results |
|---|---|---|---|
| Recognition-site scan cache (S18, `87c3fb0`, 2026-10-31) [DOCUMENTED] | 10 enzymes 100–110 s | 0.57 s | Scan-cache test |
| Cloning planner, 912 bp vector, 252 bp insert, 22 enzymes (`b4b16bb`) | 21.7 s (circular 24.0, ambiguous 23.8, no sites 4.9) | 3.2 s (3.8 / 4.3 / 2.7) | Canonical JSON of 5 planner runs plus 2 evaluations byte-identical; `test_planner_cache.py` |
| Mutation alignment memory (`f9d4356`) | 46.0 bytes/cell, 184 MB at 2,000 nt | 1.1 bytes/cell, 4.4 MB | 49 alignments identical |
| Mutation alignment time, similar sequences (`e5bb456`) | 2 kb 0.71 s; 5 kb 4.8 s; 10 kb 23 s (2 kb was 1.56 s before `f9d4356`) | 0.03 / 0.08 / 0.18 s | 300 random pairs during development plus 150 in tests identical; exactness bound |
| Sequence similarity at 10 kb (`15e2ff1`) | ~40 s under profiling (100 million cells); `compare_sequences` 10.2 s | 0.10 s | 1,200 band checks exact |
| Primer design at 10 kb (`15e2ff1`) | 135.6 s | **25.6 s** (3 kb: 9.6 s) | 3,008 hairpin results and 3 full designs identical |

**How the bands are exact.**

- **Alignment (`_banded_is_exact`).** Leaving the band |j − i| ≤ w needs at
  least G = 2w + 2 − |n − m| gap columns. Such an alignment scores at most
  match·(n + m − G)/2 + gap·G. If the banded score is strictly higher,
  every optimal path lies inside the band, and the traceback makes the
  same tie-breaking choices as the full matrix. Otherwise the band
  doubles, and the full matrix is the fallback.
- **Edit distance.** A path of cost d never leaves |i − j| ≤ d, so a banded
  result ≤ band is exact.

**Open.**

- Primer design at 10 kb is still above the audit target of 10 s.
  Profiling at 3 kb puts `_alignment_stem` (self- and heterodimer
  screening, `analyses/primer.py` around line 380) first, with 49,198
  calls and about 8 s cumulative. Next comes `normalize_dna_sequence`,
  called about 2 million times.
- Unrelated long sequences still align in O(n·m) time (about 23 s at
  10 kb).

**10 kb benchmark of all modules** (MEASURED IN SESSION, before
`15e2ff1`):

| Module | Time |
|---|---|
| Overview | 0.22 s |
| ORF | 0.07 s |
| RNA | 1.06 s |
| Protein | 0.06 s |
| Digest | 0.06 s |
| CRISPR | 0.17 s |
| Cloning planner (2 kb vector) | 3.46 s |
| Mutation | 10.2 s, since reduced to well under 1 s (similarity 0.10 s plus banded alignment) |
| Primer | 135.6 s, since reduced to 25.6 s |

**Primer3 comparison (structure screens).** The correlation is moderate
to good, but the classes overlap. A severity is a screening flag, not a
stability measure (Part III.7).

---

## Part VII – Testing, validation and quality assurance

### Current values [VERIFIED 2026-11-11]

| Check | Result |
|---|---|
| Automated tests | **421**, all passing. Run with `.venv/Scripts/python.exe -m unittest discover -s tests -t .` in about 35 s. 409 `def test_` methods plus 12 generated golden-master tests. |
| Golden master | `tests.test_golden_master` passes (12 scenarios, 13 snapshot files in `tests/golden/`). No snapshot was updated in this session. |
| Module self-tests | Run inside the suite (`test_builtin_selftests.py`). Last individual run (2026-11-03): BLAST 33/33, CRISPR 28/28, literature 64/64, mutation 44/44, variant 27/27. |
| Pyright 1.1.414 (`--pythonversion 3.10`, standard mode) | **111 errors**, no new diagnostic against the baseline of `ab4bd2a` (134) |
| Working tree | Clean |
| ReadOnly attributes | 15/15 set |

### Test families (file and method counts in Appendix D)

- **Regression.**
  - `test_p0_regressions.py`;
  - `test_phase6_tier1.py` and `test_phase6_tier2.py` (D and S items);
  - the stage tests `test_stage1`–`test_stage5`.
- **Golden master.** `tests/golden_scenarios.py` with `tests/golden/*.json`.
  Update only with `SYNAPSE_UPDATE_GOLDEN=1` after a documented scientific
  justification (Master Prompt §11.7).
- **External reference.**
  - `test_external_reference.py`: Biopython restriction (including
    Type IIS), ProtParam, IG weight;
  - `test_primer_thermodynamics.py`: primer3-py Tm;
  - `test_primer_structure_reference.py`: primer3-py ΔG.

  All are skipped when the optional package is missing.
- **Interface.**
  - `test_terminal_ui.py`: separators, `[0]`, prefixes, direct-output AST
    rules, field-message rule, layering, rendered first screens, empty
    answers;
  - `test_rendered_menus.py`: 40 scenarios, every `style.menu` call site
    (47 of 47), layout of every drawn menu, missing-file handling.
- **Architecture.** `LayeringTests` in `test_terminal_ui.py`: `core`
  imports no `ui`, also in a subprocess that refuses `ui`.
- **Workflows.** `test_integration_workflows.py` (Master Prompt 10.1–10.4
  and 10.7).
- **Smoke.** `test_main_smoke.py` runs `main.py` with scripted input.
- **New on 2026-11-11.**
  - `test_type_iis.py`;
  - `test_hgvs_notation.py`;
  - `test_rna_wobble.py`;
  - `test_orf_partial_circular.py`;
  - `test_variant_crosscheck.py`;
  - `test_planner_cache.py`;
  - `test_alignment_memory.py`;
  - `test_primer_structure_reference.py`.

### What the tests prove, and what they do not

- They prove that the implemented definitions agree with Biopython and
  Primer3 where compared.
- They prove that known defects stay fixed, that outputs are stable
  (golden), and that the menus render and react as specified.
- They do **not** prove biological correctness, and no result is
  experimentally validated.
- Golden snapshots freeze the current behavior, including any error it
  still contains.

### Gaps

- **Dialogs.** Instrumented on 2026-11-03 at `12949f3`, not re-measured:
  46 of 101 `header`, 41 of 72 `ask` and 2 of 7 `confirm` calls execute in
  the suite. The rest are export, NCBI-download, removal and
  remote-submit dialogs.
- **Coverage.** No line coverage: the `coverage` package is not installed,
  and the standard-library `trace` module is the fallback (task T7).
- **External tools.** No Python 3.10 run; no bcftools, ORFfinder,
  NEBcutter or REBASE-online comparisons; no live BLAST or literature
  runs.

### Static analysis

- Pyright runs only as `python -m pyright` (the `.exe` is blocked by
  application control).
- The baseline file is not in the repository. To reproduce it:
  1. Check out `ab4bd2a` in a temporary worktree.
  2. Set its `pyproject.toml` `venvPath` to the project folder.
  3. Run `python -m pyright --pythonversion 3.10 --outputjson`.
  4. Count the `file|rule|first message line` entries of severity
     "error". The result is 134; this was reproduced on 2026-11-11.
- The comparison ignores line numbers.

---

## Part VIII – The Phase 6 Master Prompt and the status of each stage

**File.** `C:\Users\bayra\Downloads\Project_Synapse_Phase_6_Master_Prompt_Skills_Integrated.txt`
(2,344 lines, 68,549 bytes) [VERIFIED]. It is authoritative.

**Structure.**

| § | Subject |
|---|---|
| 0 | Identity and objective |
| 1 | Access, preservation, BLAST-database policy (`data/blast_db/` deliberately excluded), external resources, skills and MCP use, approvals |
| 2 | Full audit (`docs/PHASE6_INITIAL_AUDIT.md`) |
| 3 | Scientific completion: P0 (§3.1) and per module §3.2–3.13 |
| 4 | Functional audit and FEATURE_STATUS |
| 5 | Unified terminal UI and TERMINAL_UI_STANDARD |
| 6 | Output standardization |
| 7 | Architecture |
| 8 | Language |
| 9 | STYLE_GUIDE |
| 10 | Integration workflows 10.1–10.8 |
| 11 | Testing 11.1–11.9 |
| 12 | Feature matrix |
| 13 | Stages A–I and change management |
| 14 | Boundaries (no frontend, no rewrite, no fabrication, no false claims, controlled resources) |
| 15 | Documentation |
| 16 | Final criteria |
| 17 | Final report A–K |
| 18 | Final directive |

### Stage status (assessment at HEAD, 2026-11-11)

| Stage | Objective | Done | Partial or open | Evidence |
|---|---|---|---|---|
| A Baseline | Verify copy, Git, tests; record the baseline | Done (Stage 0; initial audit `824ca5f`) | – | `docs/PHASE6_INITIAL_AUDIT.md` |
| B Critical corrections | P0 and critical defects | P0 1–9 (previous cycle); D1–D19 corrected | – | `test_p0_regressions`, tier tests, FINAL_REPORT §B |
| C P1/P2 | Resolve outstanding improvements | P1-9 (origin-spanning, 5′-partial); P2-4 (API); P2-5 (G·U count); P2-6 (planner, alignment, similarity); P2-7 (Type IIS digests) | P2-1 merge of primer generations not done; P2-2 analyses import `ui`; P2-3 CFD blocked; P2-5 "complexity" rename and RNAfold not done; P2-6 primer 10 kb at 25.6 s; P2-4 not in menus; P3 deferred | AUDIT_STATUS (not yet updated for 2026-11-11), commit messages |
| D Functional completeness | Every feature works; invalid inputs handled | Feature matrix; single-enzyme crash fixed; every menu rendered; missing files tested | Many export, download and remote dialogs not driven; restriction-site-driven primer design (10.2) absent | FEATURE_STATUS, `test_rendered_menus` |
| E Core architecture | Shared logic, API, pipeline, data models | API for every analysis; one input contract; status as data | P2-2, P2-1, pipeline menu-path parameters | LayeringTests |
| F Terminal UI | One design language | 16 interface files on `ui/style`; navigation; messages; tests | TERMINAL_UI_STANDARD §9 items | Part IX |
| G Output standardization | Results, units, exports | D18; units; 1-based positions; warnings in exports; restriction provenance | Documentation reconciliation for 2026-11-11 | Commits `317097d`..`81ccc55` |
| H Language and docs | Editorial pass, STYLE_GUIDE, README | Done on 2026-11-02 | Docs not reconciled for the 16 latest commits | `1a8f2bc`, `97086ad`, `ccf4ccb` |
| I Final verification | All tests, comparisons, final report | Performed on 2026-11-02 (`1f0c0f5`) | **Superseded**: must be redone after the remaining work; final report A–K must be rewritten | – |

Phase 6 may be declared complete only after every applicable mandatory
requirement has been checked against the Master Prompt (§16) and
documented in a new final report (§17) with the labels COMPLETE,
PARTIALLY COMPLETE, UNVERIFIED and BLOCKED.

---

## Part IX – Stage F: terminal UI standardization

### History

1. **2026-11-01.**
   - `5db21ab`: 389 separator literals in 11 widths replaced by
     `major_rule()` and `minor_rule()` at 70 columns.
   - `b1fa927`: `[0]` means only Back, Cancel or Exit; one message prefix
     per category.
   - `10431d6`: BLAST menu renumbered.
2. **2026-11-02, re-examination** (`docs/STAGE_F_UI_REVIEW.md`). It found
   that no module actually used the shared renderers.
   - Round 1 (`2e47ff9`–`623a0a8`): ORF, CRISPR, BLAST, literature, primer,
     variant.
   - Round 2 (`1ac77a4`–`f9ed98e`): session, overview, RNA, protein,
     restriction and cloning, mutation, pipeline. Also tables, report
     lines, prompts, yes/no questions, and the variant menu in §5.3 order.
3. **2026-11-03, regressions.**
   - The automatic rewrite had turned `Label [hint]:` prompts into
     `ask(default=hint)`, so Enter returned the hint and broke four
     dialogs (literature dates, primer conditions, BLAST options, cloning
     count). Fixed in `1516339`.
   - The Tm dialog lost its Enter default. Fixed in `c78800e`.
4. **2026-11-03, closing review** (`docs/STAGE_F_FINAL_REVIEW.md`) and
   recommendations:
   1. Claims corrected (`46710cd`).
   2. Direct prints converted (`4e2cace`, `e559255`, `b9ec344`):
      - 44 messages;
      - 18 headings;
      - 9 dynamic lists.
   3. FEATURE_STATUS numbers corrected (`163c54d`).
   4. Core/UI separation (`b25e800`, `b4e5a53`).
   5. Every menu rendered (`a39c39f`); this found the single-enzyme crash
      (`8eddc72`).
5. **2026-11-11.**
   - Mutation `[0] Cancel`.
   - Field-messages rule.
   - Missing-file tests.
   - Type IIS lines in the digest dialog.
   - The ORF Partial setting with three states.

### Shared renderers (`ui/style.py`, 235 lines)

| Area | Renderers |
|---|---|
| Layout | `header`, `section`, `field` / `field_line` (labels 22 wide, values at column 24), `columns`, `rule` / `major_rule` / `minor_rule` (70) |
| Menus | `menu(groups, zero=Back/Cancel/Exit/None)`, `option`, `back`, `cancel` |
| Input | `ask`, `ask_choice`, `confirm` (`[y/N]`, accepts y/yes/j/ja), `pause` |
| Messages | `error`, `warning`, `note`, `success` (ERROR:/WARNING:/NOTE:/OK:), `invalid_selection`, `notes` |

The rules are in `docs/TERMINAL_UI_STANDARD.md` §2–§9 and are enforced
by §7 items 1–11.

### Core/UI separation

- `PipelineContext.status_report()` returns the status as data.
- `print_status(render)` takes an interface's renderer, and `ui.report`
  draws it.
- A deferred (function-local) import was rejected as decoupling.

### Coverage: what "47 call sites" means

All 47 `style.menu` call sites are drawn by scenarios: 185 menu draws in
40 scenarios, plus first screens. That proves every menu renders, numbers
correctly and reacts to invalid input and 0. It does **not** cover every
dialog. Many prompts behind the menus are not driven: exports, removals,
downloads, remote submission (Part VII, gaps).

### Remaining inconsistencies (TERMINAL_UI_STANDARD §9)

- raw line input for pasted sequences;
- two CRISPR numeric helpers with their own prompts;
- closing rules;
- result lines;
- headings inside string-built reports;
- about 20 hand-padded label lines;
- the prefixes "SCIENTIFIC LIMITATION:" and "RESULT:";
- status lines such as "Loaded 2 FASTA record(s).";
- the BLAST hit selector answers "list index out of range" to an
  out-of-range number;
- the BLAST query-type question cancels on an unknown key;
- converted lists keep their own prompt labels;
- the sequence selector has no `[0]` (by design).

### P2-2

The menus still live in the analysis modules, so every analysis module
imports `ui.style`. Extracting them means:

- moving about 133 definitions (about 8,900 lines) into `ui/` modules,
  with an AST tool that collects free names;
- updating `main.py` and the tests that call the menus
  (`render_first_screen`, `test_rendered_menus`, `INTERFACE_SOURCES`, the
  direct-output rules);
- settling the exceptions: `mutation.build_mutation_report` uses
  `field_line`, `orf.select_orf` is interactive, and the self-test
  printers print.

It is a controlled, large change. Do it after T1 and T2, with
byte-identical rendered output as evidence.

---

## Part X – Development environment and safeguards

### Environment [VERIFIED 2026-11-11 unless marked]

| Item | Value |
|---|---|
| OS | Windows 11 Home (10.0.26200) |
| Shells | PowerShell 5.1 and Git Bash |
| Project copy (the only place to work) | `C:\Users\bayra\OneDrive\Desktop\Project Snyapse - Kopie`. The folder name is spelled "Snyapse". It lies in OneDrive. |
| Original project | **Never modify it.** Its path was not inspected in this session. |
| Master Prompt | `C:\Users\bayra\Downloads\Project_Synapse_Phase_6_Master_Prompt_Skills_Integrated.txt` |
| Git | 2.56.0.windows.4; user "Onur Bayraktar"; branch `main`; no remote push in this work |
| Python | 3.14.6 in `.venv` (the only interpreter; `py -0p` lists only 3.14). The project requires `>=3.10`; pyright checks at 3.10, but no 3.10 runtime exists. |
| Packages in `.venv` | biopython 1.88, primer3-py 2.4.0, numpy 2.5.4, pyright 1.1.414 (with nodejs-wheel-binaries 24.19.0, nodeenv), typing_extensions. `coverage` and `pytest` are **not** installed. |
| Runtime dependencies | None (`dependencies = []`). The optional extras are `reference = biopython, primer3-py` and `dev = pytest, pyright`. |
| External tools | BLAST+ 2.17.0 at `C:\Program Files\NCBI\blast-2.17.0+\bin`; Strawberry Perl. WSL is **not** installed, so there is no bcftools and no ORFfinder. |
| Project-local skills | `.claude/skills/`, git-ignored: biopython, bioservices, experimental-design, genomic-coordinates, gget, literature-review, primer-design, pysam, research-lookup, scientific-critical-thinking, scikit-bio |
| MCP research tools | Bio Research plugins (PubMed, bioRxiv, ChEMBL and others) appear in the session. Since 2026-11-12 they may be used under the approval policy below (no private data) |
| `.gitignore` | `__pycache__/`, `*.pyc`, `.pytest_cache/`, `/output/`, `/exports/`, `*.log`, `/data/blast_db/`, `.DS_Store`, `Thumbs.db`, `.sync/`, `/logs/`, `/.claude/skills/` |

### Line endings and file conventions [VERIFIED]

- Most source files use CRLF.
- The tests and docs use LF, as do `core/api.py` and `ui/style.py`.
- `analyses/orf.py` contains exactly two intentional bare-LF lines; keep
  them.
- Edit with exact replacements that preserve line endings. The session
  used a small helper `patchlib.apply()` that matches `\n` as `\r\n` in
  CRLF files and requires an exact match count.
- Check whitespace with `git -c core.whitespace=cr-at-eol diff --check`.

### ReadOnly protection

Exactly these 15 files carry the Windows ReadOnly attribute and must keep
it:

- `main.py`
- `analyses/`: `blast.py`, `cloning.py`, `crispr.py`, `literature.py`,
  `mutation.py`, `orf.py`, `primer.py`, `restriction.py`, `rna.py`,
  `variant_analysis.py`
- `core/pipeline.py`
- `data/`: `genetic_code.py`, `restriction_enzymes.py`
- `ui/menu.py`

To edit one, clear the attribute
(`Set-ItemProperty <file> -Name IsReadOnly -Value $false`), edit, and set
it back to `$true` at the end of the package. Count afterwards that 15 of
15 are set.

### Tooling restriction

`pyright.exe` is blocked by Windows application control. Use
`.venv\Scripts\python.exe -m pyright --pythonversion 3.10`.

### Safeguards (mandatory)

- Work only in the project copy. Never modify the original.
- Never discard pre-existing changes; check `git status` first.
- Never rewrite history, reset, delete commits or force-push.
- Preserve the ReadOnly attributes.
- Follow the approval policy below for installations, downloads and
  external requests. Never send private sequences, source code or user
  data to external services without authorization.
- Do not bypass safety restrictions or application control.
- Never change golden-master values just to make tests pass.
- Fabricate no scientific data, citations or validation.
- Preserve functionality. Prefer targeted, additive changes; never
  rewrite whole files when a targeted change suffices.
- Commit caches, secrets and skills never.
- Commit each logical package separately. End every commit message with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`, or the
  attribution the harness prescribes.

### Approval policy (owner decision of 2026-11-12)

This handover first required approval for every installation and every
external request. That was stricter than the Master Prompt. On 2026-11-12
the owner decided that the Master Prompt's rules apply (§1.3, §1.3.1,
§14.5). This supersedes the stricter wording.

**Allowed without asking**, with the purpose recorded and with the
traceability of §1.3 (source, version, date, licence, checksum):

- installing packages into the project-local `.venv` (for example
  `coverage`);
- downloading small or moderate reference data, test fixtures and
  official enzyme or sequence records, kept outside Git where large;
- consulting authoritative sources and official documentation (NCBI,
  EMBL-EBI, REBASE, PubMed, Biopython, Primer3 and others), including
  the connected MCP research tools, after checking that they work;
- live literature-provider requests with neutral queries.

**Approval required first**:

- very large downloads (for example the swissprot BLAST database, about
  0.5 GB, or complete genomes);
- system-wide installations or configuration changes (for example WSL,
  or a system-wide Python 3.10);
- anything that costs money;
- uploading private project files, source code, private sequences or
  user data. A remote BLAST of a project sequence of unknown origin
  counts as such an upload; a published reference sequence (for example
  pUC19 L09137.2) does not;
- destructive or irreversible operations.

---

## Part XI – Documentation map

| Document | Lines | Purpose | Status |
|---|---:|---|---|
| `README.md` | 181 | Purpose, modules, installation, usage, architecture, testing, limitations | Rewritten 2026-11-02. Lacks Type IIS and circular or partial ORFs. |
| `docs/CHANGES.md` | 2,028 | Change log, scientific and behavioral | Up to 2026-11-03. **Missing the 16 commits of 2026-11-11.** |
| `docs/FEATURE_STATUS.md` | 510 | Feature matrix (Master Prompt §4, §12) | Up to 2026-11-03. **Missing the 2026-11-11 changes**: cloning planner, Single-Enzyme status, Type IIS, ORF rows 6.5 and "Circular / origin-spanning ORFs" (still says not implemented), RNA hairpins, mutation notation and alignment, variant bridge, new API function, test count. |
| `docs/AUDIT_STATUS.md` | 116 | Status of every audit recommendation | Up to 2026-11-02 or 03. **P1-9, P2-4, P2-5, P2-6 and P2-7 are outdated.** |
| `docs/METHODS.md` | 565 | Algorithms, parameters, references, validation | **Current**: updated in every 2026-11-11 commit. Four references are marked as not verified at the source. |
| `docs/STYLE_GUIDE.md` | 246 | Terminology, writing, messages, code | 2026-11-02 |
| `docs/TERMINAL_UI_STANDARD.md` | 329 | UI rules and enforcement (§7 items 1–11), exceptions (§9) | Current to 2026-11-11 |
| `docs/PHASE6_INITIAL_AUDIT.md` | 940 | The audit (D/S items, evidence) | Historical (2026-10-30) |
| `docs/PHASE6_PROGRESS.md` | 164 | Numbered progress rows | Rows up to 52 (2026-11-03). **No row for 2026-11-11 yet.** |
| `docs/PHASE6_HANDOVER.md` | 488 | Handover notes of earlier sessions | Historical |
| `docs/PHASE6_FINAL_REPORT.md` | 476 | Final report A–K | Written 2026-11-02 and amended on 2026-11-03. Section J still lists items completed on 2026-11-11. To be rewritten at the end. |
| `docs/PHASE6_INDEPENDENT_REVIEW.md` | 340 | Independent review of the Phase 6 claims | Historical (2026-11-02) |
| `docs/STAGE_F_UI_REVIEW.md` | 533 | Stage F re-examination and rounds | Historical with follow-up sections 10–14 |
| `docs/STAGE_F_FINAL_REVIEW.md` | 532 | Stage F closing review and recommendations 1–5 | Historical with follow-up sections 9–10 |
| `docs/PHASE6_CONTINUATION.md` | 109 | Where the 2026-11-11 run stopped | **Current** next-step record |
| `docs/PROJECT_SYNAPSE_COMPLETE_HANDOVER.md` | this file | Full handover | Current at writing |
| `blast_results/README.md` | – | Labels the legacy BLAST result | – |

**Authority.**
1. The Master Prompt.
2. The code and tests.
3. METHODS (current).
4. PHASE6_CONTINUATION.

The remaining reports are historical snapshots.

**Inconsistencies to reconcile** (task T8):
- CHANGES, FEATURE_STATUS, AUDIT_STATUS, PHASE6_PROGRESS, FINAL_REPORT §J
  and README against the 2026-11-11 commits;
- the test count (FEATURE_STATUS says 356; the true figure is 421);
- the FINAL_REPORT pyright figure (111 is correct).

---

## Part XII – Current exact project status [VERIFIED 2026-11-11]

| Item | Value |
|---|---|
| Branch | `main` |
| HEAD at the start of this documentation task | `690c9e8` "Correct the test count in the continuation record (421)". The commit that adds this file follows it; check with `git log --oneline -3`. Neither `12949f3` nor `690c9e8` is the HEAD after this file is committed. |
| Commits | 125 before this file (Appendix A) |
| Working tree | Clean before this file |
| Tests | 421 run, OK |
| Golden master | OK |
| Pyright | 111 errors, no new diagnostic against the 134 baseline |
| ReadOnly | 15/15 |
| Latest continuation record | `docs/PHASE6_CONTINUATION.md` (commits `b8e92d8`, `690c9e8`) |
| Latest progress row | `docs/PHASE6_PROGRESS.md` row 52 (2026-11-03) |

### Package A – test gaps and navigation (2026-11-11)

| Commit | Content | Validation |
|---|---|---|
| `68e59d1` | HTTPError responses closed in `literature._request_text` and `NCBIBlastClient._request` | Two tests keep the errors alive; both fail on the old code; no ResourceWarning remains |
| `8379ac2` | Mutation `[0] Cancel` via `SequencePromptCancelled`; 10 field-messages converted to ERROR/NOTE lines; new rule `test_no_message_is_drawn_as_a_field`; `MissingFileTests` (loaders and five dialogs) | The rule finds all 10 sites in the old code; rendered scenario `mutation_cancel` |

### Package B – science and performance (2026-11-11)

| Commit | Content |
|---|---|
| `b4b16bb` | Cloning planner search cache. The S18 test now expects 0 repeated scan requests. |
| `4939688` | Type IIS enzymes; digest dialog shows them |
| `b7eb1de` | HGVS `p.(...)` notation |
| `1ae2833` | Restriction coordinate convention |
| `fedb750` | G·U wobble count |
| `f020a1d` | Origin-spanning and 5′-partial ORFs |
| `27255df` | primer3-py structure comparison |
| `2c68e61` | P2-4 variant cross-check |
| `f9d4356` | One-byte traceback in the alignment |
| `e5bb456` | Exact banded alignment |
| `15e2ff1` | Banded Levenshtein and faster hairpin screen |

Then the stop record: `b8e92d8`, `690c9e8`.

### The next task (the plan is still current)

The order in `docs/PHASE6_CONTINUATION.md` was checked against the code:

1. Primer design under 10 s at 10 kb (`_alignment_stem`).
2. Primer design driven by a chosen restriction site.
3. P2-2.
4. Reduce pyright errors.
5. Coverage with the standard-library `trace` module.
6. **Tool and data ledger.** This is a documentation table of versions,
   sources, dates, licences and checksums (Master Prompt 1.3.1, 11.8), not
   the creation of directories.
7. Reconcile the documentation for all 16 commits of 2026-11-11.
8. Full verification, independent audit and final report A–K.

The continuation record also lists P2-1 (merging the primer generations;
decide whether it is required) and the pipeline menu-path parameters for
sequence overview and cloning.

---

## Part XIII – Task register, blockers and approvals

| ID | Priority | Task | Why open | Files | Approach | Verification | Approval |
|---|---|---|---|---|---|---|---|
| T1 | 1 | Primer design 10 kb < 10 s (P2-6) | 25.6 s now | `analyses/primer.py` (`_alignment_stem` ~l.380; `_make_candidate` ~l.509; `normalize_dna_sequence` ~l.223) | 1. Snapshot `design_primers_professional` on 600/900/1,200 bp (serialized dataclasses) and the structure results of about 3,000 primers. 2. Optimize the dimer stem search (precompute reverse complements, avoid re-normalization, early exits that cannot change results). 3. Prove the output identical. | Snapshots identical; golden; 10 kb timing; suite | No |
| T2 | 2 | Restriction-site-driven primer design (Master Prompt 10.2) | No wiring exists | `analyses/primer.py`, `analyses/restriction.py`, `core/api.py`, primer and restriction menus | A chosen site (enzyme and position) defines the target or flank region passed to `design_primers`; record it in the parameters | Workflow test 10.2; rendered dialog | No |
| T3 | 3 | P2-2: menus out of the analysis modules | Analyses import `ui.style` | all 12 analysis modules, `ui/`, `main.py`, tests | AST-guided move of style-using definitions to `ui/<module>_menu.py`; update callers; test that imports every analysis module with `ui` refused | Byte-identical rendered output (`test_rendered_menus`), suite, golden, layering test | No (large; analyze first) |
| T4 | 3 | Pipeline menu-path parameters (sequence overview, cloning) | Which inputs to record is not yet defined | `core/pipeline.py`, `main.py` | Name the parameters explicitly in `store()` | Report export test | No |
| T5 | 4 | P2-1: merge the surviving primer generations | Deliberately deferred | `analyses/primer.py` | Decide first; if done, prove with bytecode and golden checks | Golden, snapshots | No |
| T6 | 4 | Pyright triage (111) | Reachable `None` paths | primer, CRISPR, mutation, tests | Fix real defects first; no blanket ignores | `pyright_diff`: fewer errors, nothing new | No |
| T7 | 4 | Line coverage | `coverage` not installed | – | Standard-library `trace` over the suite; report the uncovered dialogs | Report in docs | No (since 2026-11-12 a project-local `coverage` install is allowed) |
| T8 | 5 | Documentation reconciliation for 2026-11-11 | Not yet done | CHANGES, FEATURE_STATUS, AUDIT_STATUS, PHASE6_PROGRESS, FINAL_REPORT §J, README | One entry or row per commit | Table checker; cross-check against the code | No |
| T9 | 5 | Tool and data ledger | Versions recorded, but not sources, dates, licences or checksums | new section in METHODS or a ledger doc | List Biopython 1.88, primer3-py 2.4.0, BLAST+ 2.17.0, pyright 1.1.414, pUC19 L09137.2 and the fixtures with SHA-256 | – | No |
| T10 | 6 | Full verification, independent audit, final report A–K | Final stage | all docs | Re-read the Master Prompt; build the requirement matrix COMPLETE / PARTIAL / UNVERIFIED / BLOCKED | – | No |
| T11 | – | Remaining Stage F §9 items (hand-padded labels, BLAST hit selector reply, prefixes) | Minor | various | Targeted fixes with tests | Rendered tests | No |
| T12 | – | P2-4 wired into the menus | No record of the VCF contig mapped to the sequence | variant UI, pipeline | Store the contig mapping when normalizing, then offer the check in primer and CRISPR | Tests | No |
| T13 | – | Banded alignment for unrelated long sequences | O(n·m) remains | `analyses/mutation.py` | Probably acceptable; document | – | No |

### Blocked (need the owner's decision)

Re-assessed on 2026-11-12 under the approval policy of Part X.

| Item | Needs |
|---|---|
| Python 3.10 runtime test | A 3.10 interpreter. A system-wide install needs approval; a project-local interpreter kept outside Git may be used if one can be obtained from an official source |
| `bcftools norm` and ORFfinder comparisons | WSL or a Linux host (system-wide: approval) |
| Live NCBI BLAST run and RID resume | Allowed with a published reference sequence; a project sequence of unknown origin needs approval |
| Live literature providers (S11) | Allowed with neutral queries since 2026-11-12 |
| REBASE online re-check, NEBcutter, PubMed check of four METHODS references, the Guruprasad 1990 primary table | Allowed since 2026-11-12 (consulting public sources); record source and date |
| CFD off-target scoring (P2-3) | The published matrix as a data file with citation and licence; obtaining it from the publication is allowed |
| `coverage` package | Allowed since 2026-11-12 (project-local) |
| Licence (`pyproject.toml`: "decision pending"); `blast_results/TEST_8_REFERENCE_blast.json` | Owner decisions |
| Local BLAST on a real database (swissprot about 0.5 GB) | Approval for the download |
| P3-1 Experiment planner, P3-2 genome specificity, P3-3 front end | Deferred by design; P3-3 is forbidden in Phase 6 |

No additional blockers were found during this inspection.

---

## Part XIV – Working with the project owner

- **Language.** Speak German with the owner. Write code, comments,
  commit messages and technical documentation in English.
- **Honesty.** Be factual; never overstate quality or completion; say what
  was not verified.
- Explain scientific concepts when they matter. Give copy-paste-ready
  commands or code when asked.
- **Autonomy.** Work in logical packages once authorized, and do not ask
  for confirmation after every small step. Still ask before anything that
  needs approval (Part X).
- **Continuity.** Do not restart the Master Prompt after a context change;
  continue from the verified HEAD and the continuation record.
- Preserve architectural decisions and existing functionality. Explain
  why each change is needed. Test before claiming success. Be especially
  careful with scientific correctness.
- **Recurring owner expectations seen in this conversation:**
  - independent reviews of claims;
  - tests that fail on the old code;
  - byte-identical evidence for behavior-preserving changes;
  - stating which items are BLOCKED and why;
  - stopping and reporting when asked;
  - keeping the context budget in mind; the owner asked to stop after
    package B for that reason.
- **Goal.** A serious, scientifically credible bioinformatics platform, not
  a demonstration.

---

## Part XV – Exact instructions for the next session

1. Read this document completely.
2. Read the Master Prompt (path in Part X) completely. It is authoritative.
3. Read `docs/PHASE6_CONTINUATION.md`.
4. Run `git log --oneline -5` and `git status`. Do not assume the HEAD or
   test counts given here; re-verify with the commands in Appendix E,
   including ReadOnly 15/15 and pyright against the 134 baseline. The
   baseline file must be regenerated (Part VII).
5. Verify the project state. If it differs from Part XII, find out why
   before changing anything.
6. Continue with task **T1**, then T2 to T10 in the order of Part XIII,
   unless the owner says otherwise.
7. Do not repeat completed work (Part XII, packages A and B).
8. Keep every safeguard (Part X).
9. Work in logical packages. After each package:
   - run the relevant and the full tests, golden master and pyright;
   - inspect the diff and the line endings;
   - restore ReadOnly;
   - commit with the attribution trailer;
   - update `docs/PHASE6_CONTINUATION.md`.
10. Phase 6 is **not complete**. Declare it complete only after T10, with
    every Master Prompt requirement checked and labeled.

Do not start unrelated features or any frontend or UI redesign. Do not
restart Phase 6.

---

## Part XVI – Completeness check of this document

Checked on 2026-11-11 against:

- the complete Master Prompt (read in full in this session);
- `git log` (125 commits; Appendix A);
- the source tree and file sizes (Appendix B);
- `core/api.py` `__all__` (Appendix C);
- the test files (Appendix D);
- README, METHODS, TERMINAL_UI_STANDARD, the continuation record,
  FINAL_REPORT §B and §J, AUDIT_STATUS P-items, PHASE6_PROGRESS and the
  `.gitignore`.

Fresh runs: suite, golden master, pyright, ReadOnly count, `git status`.

**Not independently verifiable:**
- the content and quality of the previous "Claude Cowork" cycle beyond its
  commits;
- the original project's location;
- the owner's motivation beyond the Master Prompt's own text.

**Not re-measured:**
- all benchmark timings (Part VI);
- the 2026-11-03 dialog-coverage figures;
- the individual self-test counts.

**Contradictions resolved:**
- the plan item "create tool and data directories" means the tool and
  data ledger;
- HEAD is newer than both `12949f3` and `690c9e8`.

---

## Appendix A – Complete commit history [VERIFIED 2026-11-11]

All 125 commits, oldest first. The commit that adds this file follows them.

| # | Commit | Date | Subject |
|---:|---|---|---|
| 1 | `e4395ed` | 2026-10-29 | Baseline: project state before audit-driven changes |
| 2 | `f4a2ac0` | 2026-10-29 | Stage 0: automated test suite, golden master, P0 reproduction tests |
| 3 | `08fdce2` | 2026-10-29 | Fix Python 3.10/3.11 import error in protein.py |
| 4 | `59d1997` | 2026-10-29 | Stage 1: fix the nine P0 audit findings |
| 5 | `b99e63d` | 2026-10-29 | Stage 2: shared data model, visible result storage, file input, logging |
| 6 | `5171926` | 2026-10-29 | Stage 3: P1 scientific improvements and reproducible report |
| 7 | `dd15af5` | 2026-10-29 | Stage 4: consolidation, shared tables, core API, methods documentation |
| 8 | `88cc4ae` | 2026-10-29 | Stage 5: unresolved residues keep their position in protein analysis |
| 9 | `0dad3bb` | 2026-10-29 | Add docs/AUDIT_STATUS.md: assessment of every audit recommendation |
| 10 | `5242fcb` | 2026-10-30 | Record pyright type-checking configuration from the previous cycle |
| 11 | `5d8c73a` | 2026-10-30 | Ignore locally installed Claude skills (/.claude/skills/) |
| 12 | `824ca5f` | 2026-10-30 | Add Phase 6 initial audit and feature verification matrix |
| 13 | `d4544ad` | 2026-10-31 | Flag variants whose contig has no reference instead of skipping them (D12) |
| 14 | `5c098e3` | 2026-10-31 | Use one restriction-site matcher; ambiguity never makes a cloning site (D4) |
| 15 | `47146b0` | 2026-10-31 | Honor circular topology in the restriction map and site results (D3) |
| 16 | `7334296` | 2026-10-31 | Honor vector and insert topology in cloning analysis (D5) |
| 17 | `dc780f6` | 2026-10-31 | Decode initiator Met only for initiation codons of the genetic code (D1) |
| 18 | `ab4bd2a` | 2026-10-31 | Translate an ORF with its own genetic code (D2) |
| 19 | `02d2899` | 2026-10-31 | Add Phase 6 progress record and the D14 sub-findings |
| 20 | `b8f5ea8` | 2026-10-31 | Use one BLAST query-type vocabulary (D14-F2) |
| 21 | `d057ac1` | 2026-10-31 | Add central sequence alphabets and a type classifier (D14, D13) |
| 22 | `776de4a` | 2026-10-31 | Read, keep and report the declared molecule type (D14-F3, F11) |
| 23 | `7c716a2` | 2026-10-31 | Take the BLAST query type from the central classifier (D14-F1, S20) |
| 24 | `162ca16` | 2026-10-31 | One input contract for menu and API; route protein records (D13, D14) |
| 25 | `343ffb7` | 2026-10-31 | Keep IUPAC protein codes in place and give them meaning (D14-F6) |
| 26 | `d40d328` | 2026-10-31 | Normalize and check input in ORF search, transcription, GC count (D14) |
| 27 | `5b120ec` | 2026-10-31 | Address a local BLAST database by its path prefix (D6, S7) |
| 28 | `e1a48c6` | 2026-10-31 | Label the stale BLAST result file as a legacy record (D7) |
| 29 | `7a4a718` | 2026-10-31 | Add handover notes to the Phase 6 progress record |
| 30 | `3514617` | 2026-10-31 | Do not count a guide's own site as an off-target (D9) |
| 31 | `ca2ec39` | 2026-10-31 | Report the modeled CRISPR cut site with its source (D8) (subject written in American spelling here; the commit itself uses the British double-l form) |
| 32 | `8145556` | 2026-10-31 | Make the CRISPR config fields act; keep exports in the project (D10, D11) |
| 33 | `e47cb85` | 2026-10-31 | Let the pipeline's ORF runner take the menu's search settings (D19) |
| 34 | `f3a09c1` | 2026-10-31 | Label the unimplemented AI Experiment Planner as such (D15) |
| 35 | `c71963b` | 2026-10-31 | Fix the variant-to-mutation bridge and a flag leak (S15, S16) |
| 36 | `c743beb` | 2026-10-31 | Cut the cloning digest at the enzyme's cut, not its site start (S17) |
| 37 | `6c20e69` | 2026-10-31 | Classify and score ORF start codons with the genetic code (S19) |
| 38 | `642f31b` | 2026-10-31 | Add Phase 6 handover report for the planned pause |
| 39 | `6f7987b` | 2026-10-31 | Correct the guide seed end and the RTT coverage check (S23, S3, S4, S5) |
| 40 | `bff3c7d` | 2026-10-31 | Record the restricted CRISPR assignment in the handover report |
| 41 | `c52f9bb` | 2026-10-31 | Give the primer Tm model one name instead of shadowing it (D16) |
| 42 | `5a3220e` | 2026-10-31 | Report what the GenBank reader, the serializer and the REF check discard (S8, S10, S21, S22) |
| 43 | `8c47b17` | 2026-10-31 | Assess five suspected findings; quantify the alignment cost (S1, S2, S9, S11, S12) |
| 44 | `edbc28f` | 2026-10-31 | Update the handover report for the second Opus 5 assignment |
| 45 | `dfe773d` | 2026-10-31 | Remove the superseded primer generations after a reachability review (D16 part 2) |
| 46 | `8d28586` | 2026-10-31 | Use the selected ORF's genetic code in the RNA codon view (P1-9) |
| 47 | `87c3fb0` | 2026-10-31 | Cache recognition-site scans so the planner stops rescanning (S18) |
| 48 | `6a8c00c` | 2026-10-31 | Give the API entry points for the remaining analyses (S14) |
| 49 | `ad843c5` | 2026-10-31 | Resume an NCBI BLAST search from its request ID (S6, P1-7) |
| 50 | `7008cd8` | 2026-10-31 | Run the session from main() instead of at import (S13) |
| 51 | `692cb20` | 2026-10-31 | Update the handover report for the third Opus 5 assignment |
| 52 | `5db21ab` | 2026-11-01 | Draw every terminal separator at one width through ui.style (Stage F, D17) |
| 53 | `b1fa927` | 2026-11-01 | Give [0] one meaning and messages one prefix; document and enforce the UI standard (Stage F, D17) |
| 54 | `10431d6` | 2026-11-01 | Renumber the BLAST menu so its groups read in order (Stage F) |
| 55 | `f803577` | 2026-11-01 | Measure the XML parsing exposure and classify every outstanding item (S24) |
| 56 | `aa5703b` | 2026-11-01 | Update the handover report for Stage F |
| 57 | `317097d` | 2026-11-02 | Stop calling heuristic scores "confidence" in data and exports (Stage G, D18) |
| 58 | `bc89394` | 2026-11-02 | Report the sequence-type outcome as an evidence level (Stage G, D18) |
| 59 | `ebf4d52` | 2026-11-02 | Use one notation per unit and show primer positions 1-based (Stage G) |
| 60 | `cef464a` | 2026-11-02 | Carry a result's own warnings into the exported summary (Stage G, 6.4) |
| 61 | `81ccc55` | 2026-11-02 | Record restriction parameters and the enzyme data version (Stage G, 6.4) |
| 62 | `1a8f2bc` | 2026-11-02 | Edit text for one spelling, one wording per message and a sober tone (Stage H) |
| 63 | `4e5dbdc` | 2026-11-02 | Print every BLAST E-value with three significant figures (Stage G/H) |
| 64 | `97086ad` | 2026-11-02 | Add the style guide and check its enforcement table (Stage H) |
| 65 | `ccf4ccb` | 2026-11-02 | Rewrite the README to cover what section 15.1 asks for (Stage H) |
| 66 | `1269335` | 2026-11-02 | State what "sequence complexity" measures and where it fails (Stage H, P2-5) |
| 67 | `4c31602` | 2026-11-02 | Record the Stage G and H checkpoint in the handover |
| 68 | `2851e11` | 2026-11-02 | Compare primer Tm with primer3-py and record the result (Stage I, P1-3) |
| 69 | `835744e` | 2026-11-02 | Correct the IG instability weight found by new Biopython reference tests (Stage I) |
| 70 | `999f77a` | 2026-11-02 | Test the cross-module workflows end to end (Stage I, 10.1-10.4, 10.7) |
| 71 | `82ae2e9` | 2026-11-02 | Bring AUDIT_STATUS up to date for the end of Phase 6 (Stage I, 15.4) |
| 72 | `70cf362` | 2026-11-02 | Bring the feature matrix to the end-of-Phase-6 state (Stage I) |
| 73 | `db74504` | 2026-11-02 | Make the CRISPR module importable on Python 3.10 again (Stage I, 11.8) |
| 74 | `1f0c0f5` | 2026-11-02 | Add the Phase 6 final report and close the handover (Stage I) |
| 75 | `fc77ff4` | 2026-11-02 | Add the Phase 6 closing review |
| 76 | `fc13aad` | 2026-11-02 | Correct the Phase 6 documents where the closing review found errors |
| 77 | `2e47ff9` | 2026-11-02 | Route the ORF screens through ui.style (Stage F completion, part 1) |
| 78 | `e7eaf58` | 2026-11-02 | Route the CRISPR screens through ui.style (Stage F completion, part 2) |
| 79 | `5cb4c9a` | 2026-11-02 | Route the BLAST screens through ui.style (Stage F completion, part 3) |
| 80 | `b181ae3` | 2026-11-02 | Route the literature screens through ui.style (Stage F completion, part 4) |
| 81 | `fb4070a` | 2026-11-02 | Route the primer screens through ui.style (Stage F completion, part 5) |
| 82 | `fdc968d` | 2026-11-02 | Route the variant screens through ui.style (Stage F completion, part 6) |
| 83 | `b80fd37` | 2026-11-02 | Wait for Enter after CRISPR precondition messages (Stage F completion) |
| 84 | `623a0a8` | 2026-11-02 | Test the rendered screens and record the Stage F completion |
| 85 | `1ac77a4` | 2026-11-02 | Route the session screens through ui.style (Stage F completion, part 7) |
| 86 | `bf6ae87` | 2026-11-02 | Route the sequence overview through ui.style (Stage F completion, part 8) |
| 87 | `9b7f2ef` | 2026-11-02 | Route the RNA and protein screens through ui.style (Stage F completion, part 9) |
| 88 | `85c22ac` | 2026-11-02 | Route the restriction and cloning screens through ui.style (Stage F completion, part 10) |
| 89 | `a4df3d4` | 2026-11-02 | Route the mutation screens and pipeline status through ui.style (Stage F completion, part 11) |
| 90 | `3853d87` | 2026-11-02 | One form for result tables, report lines, prompts and yes/no questions (Stage F completion, part 12) |
| 91 | `69e3f0f` | 2026-11-02 | Order the variant menu as the work is done (Stage F completion, part 13) |
| 92 | `7970e5a` | 2026-11-02 | Test every interface screen as rendered (Stage F completion, part 14) |
| 93 | `7ecb88d` | 2026-11-02 | Draw the ORF selection and topology lists through style.menu (Stage F completion, part 15) |
| 94 | `f9ed98e` | 2026-11-02 | Record the second Stage F round in the standard, the review and the reports |
| 95 | `1516339` | 2026-11-03 | Let Enter keep the default in six dialogs again (Stage F regression) |
| 96 | `46710cd` | 2026-11-03 | Record the Stage F closing review and correct the overstated UI claims |
| 97 | `4e2cace` | 2026-11-03 | Print the remaining validation and status messages through ui/style (Stage F, recommendation 2a) |
| 98 | `e559255` | 2026-11-03 | Draw the remaining directly printed headings through style.section (Stage F, recommendation 2b) |
| 99 | `b9ec344` | 2026-11-03 | Draw the dynamic BLAST and CRISPR option lists through style.menu (Stage F, recommendation 2c) |
| 100 | `73a1d93` | 2026-11-03 | Detect computed menu keys, plain messages and plain headings in every interface file |
| 101 | `163c54d` | 2026-11-03 | Bring the menu numbers and summary figures in FEATURE_STATUS in line with the code |
| 102 | `ccc987b` | 2026-11-03 | Type the f-string walk in the direct-output rule so pyright reports nothing new |
| 103 | `228b800` | 2026-11-03 | Record recommendations 1-3 of the Stage F closing review in the standard and the reports |
| 104 | `b25e800` | 2026-11-03 | Let core hand the pipeline status to a renderer instead of importing ui (Stage F, recommendation 4) |
| 105 | `b4e5a53` | 2026-11-03 | Test that core imports and runs with the ui package unavailable |
| 106 | `c78800e` | 2026-11-03 | Let Enter choose the nearest-neighbor Tm model again (Stage F regression) |
| 107 | `61d0d38` | 2026-11-03 | Print three more input-check messages as ERROR lines |
| 108 | `8eddc72` | 2026-11-03 | Give the single-enzyme cloning analysis a report of its own (it always crashed) |
| 109 | `a39c39f` | 2026-11-03 | Render and drive every menu of the interface (Stage F, recommendation 5) |
| 110 | `12949f3` | 2026-11-03 | Record recommendations 4 and 5 in the standard, the reviews and the reports |
| 111 | `68e59d1` | 2026-11-11 | Close the response an HTTPError carries in the literature and BLAST clients |
| 112 | `8379ac2` | 2026-11-11 | Give the mutation sequence prompt a cancel, report errors as ERROR lines, test missing files |
| 113 | `b4b16bb` | 2026-11-11 | Compute each enzyme analysis once per cloning-planner search (P2-6) |
| 114 | `4939688` | 2026-11-11 | Digest with Type IIS enzymes on both strands (P2-7) |
| 115 | `b7eb1de` | 2026-11-11 | Write predicted protein changes as HGVS p.(...) and use them in reports and variants |
| 116 | `1ae2833` | 2026-11-11 | Use one coordinate convention for restriction site and cut records |
| 117 | `fedb750` | 2026-11-11 | Count G·U wobble pairs in the RNA hairpin screen as a separate figure (P2-5) |
| 118 | `f020a1d` | 2026-11-11 | Report ORFs across the origin of circular sequences and 5'-partial ORFs (P1-9) |
| 119 | `27255df` | 2026-11-11 | Compare the primer hairpin and dimer screens with primer3-py |
| 120 | `2c68e61` | 2026-11-11 | Check designed primers and guides against variants (P2-4) |
| 121 | `f9d4356` | 2026-11-11 | Keep the mutation alignment to one traceback byte per cell (P2-6) |
| 122 | `e5bb456` | 2026-11-11 | Align large similar sequences in an exact band (P2-6, 10 kb under 10 s) |
| 123 | `15e2ff1` | 2026-11-11 | Speed up sequence similarity and the primer hairpin screen without changing results (P2-6) |
| 124 | `b8e92d8` | 2026-11-11 | Record where the Phase 6 completion run stopped and what remains |
| 125 | `690c9e8` | 2026-11-11 | Correct the test count in the continuation record (421) |

## Appendix B – Source files and sizes [VERIFIED 2026-11-11]

| Lines | File |
|---:|---|
| 0 | `analyses/__init__.py` |
| 0 | `core/__init__.py` |
| 0 | `data/__init__.py` |
| 0 | `fileio/__init__.py` |
| 0 | `ui/__init__.py` |
| 4 | `data/genetic_code.py` |
| 26 | `core/validation.py` |
| 37 | `core/logging_setup.py` |
| 63 | `ui/menu.py` |
| 78 | `fileio/fasta.py` |
| 107 | `fileio/loader.py` |
| 121 | `data/amino_acids.py` |
| 146 | `core/orf_context.py` |
| 172 | `core/seqtools.py` |
| 204 | `core/report.py` |
| 207 | `fileio/genbank.py` |
| 233 | `core/models.py` |
| 235 | `ui/style.py` |
| 295 | `core/alphabets.py` |
| 421 | `data/restriction_enzymes.py` |
| 562 | `core/api.py` |
| 721 | `ui/report.py` |
| 938 | `main.py` |
| 1317 | `core/pipeline.py` |
| 2395 | `analyses/literature.py` |
| 2493 | `analyses/sequence.py` |
| 2830 | `analyses/restriction.py` |
| 2910 | `analyses/variant_analysis.py` |
| 3083 | `analyses/crispr.py` |
| 3093 | `analyses/rna.py` |
| 3360 | `analyses/orf.py` |
| 4369 | `analyses/protein.py` |
| 5057 | `analyses/cloning.py` |
| 5058 | `analyses/primer.py` |
| 5273 | `analyses/mutation.py` |
| 5576 | `analyses/blast.py` |

Total: 51384 lines of Python in the application packages. Other folders: `examples/` (dna.fasta, crispr_batch_test.fasta), `tests/fixtures/` (standard_code_reference.py, synthetic_circular.gb), `tests/golden/` (13 JSON snapshots), `LEGACY/` (first_dna_analyzer.py, dna.txt, a screenshot), `Synapse_Website/` (out of scope), `blast_results/` (legacy record), `exports/` and `logs/` (git-ignored).

## Appendix C – Public API (`core/api.py`, `__all__`) [VERIFIED]

`load_records`, `sequence_overview`, `find_orfs`, `analyze_rna`, `analyze_protein_from_cds`, `analyze_protein`, `restriction_digest`, `design_primers`, `compare_sequences`, `analyze_vcf`, `design_crispr_guides`, `plan_cloning`, `filter_variants`, `mutation_report`, `blast_search`, `blast_resume`, `literature_search`, `cross_check_variants` (since 2026-11-11). Each returns an `AnalysisResult` (module, data, parameters, warnings, limitations, input hash) except `cross_check_variants`, which returns a list of findings. `_sequence_of()` enforces the input contract (type check, declared type, topology).

## Appendix D – Test files [VERIFIED 2026-11-11]

| `def test_` methods | File |
|---:|---|
| 90 | `test_phase6_tier2.py` |
| 46 | `test_phase6_sequence_types.py` |
| 42 | `test_phase6_tier1.py` |
| 31 | `test_terminal_ui.py` |
| 28 | `test_rendered_menus.py` |
| 23 | `test_stage3_additions.py` |
| 22 | `test_p0_regressions.py` |
| 15 | `test_stage1_additions.py` |
| 15 | `test_orf_partial_circular.py` |
| 14 | `test_stage2_foundation.py` |
| 11 | `test_stage4_consolidation.py` |
| 10 | `test_type_iis.py` |
| 9 | `test_variant_crosscheck.py` |
| 9 | `test_alignment_memory.py` |
| 8 | `test_primer_thermodynamics.py` |
| 8 | `test_external_reference.py` |
| 7 | `test_integration_workflows.py` |
| 5 | `test_stage5_protein_positions.py` |
| 4 | `test_hgvs_notation.py` |
| 3 | `test_rna_wobble.py` |
| 3 | `test_primer_structure_reference.py` |
| 3 | `test_planner_cache.py` |
| 2 | `test_main_smoke.py` |
| 1 | `test_builtin_selftests.py` |
| 0 | `test_golden_master.py` |

`test_golden_master.py` generates 12 tests from `tests/golden_scenarios.py`, so the total is 409 + 12 = 421. Support files: `tests/_support.py` (`deterministic_dna`, `reverse_complement`, `PROJECT_ROOT`, golden helpers), `tests/fixtures/`.

## Appendix E – Commands

Run from the project root.

**Git Bash**

```bash
.venv/Scripts/python.exe -m unittest discover -s tests -t .            # full suite (~35 s)
.venv/Scripts/python.exe -m unittest tests.test_golden_master          # golden master
SYNAPSE_UPDATE_GOLDEN=1 .venv/Scripts/python.exe -m unittest tests.test_golden_master   # ONLY after a justified, documented change
.venv/Scripts/python.exe -m unittest tests.test_rendered_menus tests.test_terminal_ui
.venv/Scripts/python.exe -m pyright --pythonversion 3.10               # never pyright.exe
.venv/Scripts/python.exe -m core.pipeline                              # core smoke test, prints JSON
.venv/Scripts/python.exe main.py                                       # interactive session
git -c core.whitespace=cr-at-eol diff --check
```

**PowerShell: check ReadOnly (expect 15)**

```powershell
$f='main.py','analyses\blast.py','analyses\cloning.py','analyses\crispr.py','analyses\literature.py','analyses\mutation.py','analyses\orf.py','analyses\primer.py','analyses\restriction.py','analyses\rna.py','analyses\variant_analysis.py','core\pipeline.py','data\genetic_code.py','data\restriction_enzymes.py','ui\menu.py'
($f | Where-Object { (Get-Item $_).IsReadOnly }).Count
```

**Rebuild the pyright baseline (134)**

1. `git worktree add --detach <tmp> ab4bd2a`.
2. In `<tmp>/pyproject.toml`, set `venvPath` to the project folder.
3. Run `python -m pyright --pythonversion 3.10 --outputjson` there.
4. Write the `file|rule|first message line` entries of severity "error"
   to a file.
5. `git worktree remove --force <tmp>`.

## Appendix F – Status vocabulary and glossary

- **Feature statuses (Master Prompt §4.1).** Implemented, Execution-tested,
  Integration-tested, Independently validated, Partially implemented,
  Not implemented, Blocked by external dependency, Not applicable.
  "Reachable" in FEATURE_STATUS means read but not driven.
- **Final-report labels (owner's instruction).** COMPLETE, PARTIALLY
  COMPLETE, UNVERIFIED, BLOCKED. Never write "complete" while a mandatory
  item is open.
- **Item identifiers.**
  - **P0-x, P1-x, P2-x, P3-x:** recommendations of the previous audit (see
    AUDIT_STATUS).
  - **D1–D19:** confirmed Phase 6 defects.
  - **S1–S24:** suspected items.
  - **F-items (D14-F1 … F12):** sub-findings of D14.
- **Site status.**
  - **CONFIRMED:** a site of exact bases.
  - **POSSIBLE:** a site that matches only through IUPAC ambiguity.
- **Golden master.** Frozen JSON outputs of fixed scenarios; a change
  needs a scientific justification (§11.7).
- **Pipeline.** `SynapsePipeline` and `PipelineContext`: the session's
  result store.
- **StatusRenderer.** A callable that receives
  `PipelineContext.status_report()`.
- **Band, exactness bound.** See Part VI.
- **Wobble pair.** A G·U pair in RNA.
- **5′-partial and 3′-partial ORF.** An open frame without a start codon in
  the sequence, or without a stop codon.
- **Origin-spanning ORF.** An ORF on a circular sequence whose
  coordinates run past the origin; it carries `wraps_origin` and
  `stop_on_circle`.
