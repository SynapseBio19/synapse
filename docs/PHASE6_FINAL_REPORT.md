# Phase 6 final report

Project Synapse, Phase 6, 2026-10-30 to 2026-11-13. Working copy
`Project Snyapse - Kopie`; the original project was not modified.
Governing specification: the Phase 6 Master Prompt (authoritative over
every report). Work by Claude Opus 5 and Claude Opus 5.5.

This report replaces the version of 2026-11-02/03 (in the Git history at
`12949f3`), which a later review found to overstate completion and which
the completion run of 2026-11-11 to 2026-11-13 superseded. Status labels:
**COMPLETE**, **PARTIALLY COMPLETE**, **UNVERIFIED**, **BLOCKED**. The word
COMPLETE is used for a requirement only when everything it demands is done
and checked.

Detailed records: `docs/PHASE6_INITIAL_AUDIT.md` (findings D1-D19, S1-S24),
`docs/PHASE6_PROGRESS.md` (one row per work package), `docs/CHANGES.md`
(every change with before, after and reason), `docs/FEATURE_STATUS.md`
(every feature), `docs/AUDIT_STATUS.md` (every audit recommendation),
`docs/METHODS.md` (methods and validation), `docs/TEST_COVERAGE.md`,
`docs/TOOL_AND_DATA_LEDGER.md`, `docs/PHASE6_CONTINUATION.md`.

---

## A. Executive summary

**Outcome.** Synapse is a coherent terminal application whose science no
longer depends on its interface. Every analysis runs through `core.api`
with `ui` unimportable; the eleven module dialogs live in `ui/` and share
one renderer. The defects confirmed during Phase 6 (D1-D25, the IG
instability weight and the Stage F faults: the single-enzyme crash, two
Enter regressions, the BLAST hit selector and query-type question) are
corrected, each with a regression test. The Master Prompt's open scientific items
(P1-9, P2-2, P2-4, P2-6, P2-7) are done; the 10 kb performance target is
met for every module.

**Readiness.** Research and teaching software. Its calculations agree with
established implementations where compared (Biopython, primer3-py,
REBASE, NCBI BLAST); none of its results is experimentally validated.

**Main improvements of the final run (2026-11-11 to 2026-11-13).**
- Primer design at 10 kb: 135.6 s -> 9.6 s CPU time, identical results.
- Primers designed around a chosen restriction site (workflow 10.2).
- Menus out of the analysis modules (P2-2), 55 screens byte-identical.
- Type IIS enzymes; circular and 5'-partial ORFs; HGVS-style notation;
  variant cross-check of primers and guides in the menus.
- Six further defects (D20-D25), among them a structural penalty that
  never applied, a primer risk assessment missing from every export, and
  a network error that lost a submitted BLAST search.
- Suite on Python 3.10.11 and 3.14.6; coverage measured (71 %); every
  enzyme definition checked against REBASE itself; live checks of BLAST
  and the literature providers.

**Verdict.** Phase 6 is **not declared COMPLETE**. All safely executable
work is done and verified. Open are: 10 items **BLOCKED** (they need the
owner's approval, system-wide tools, very large data, a licensed data
file or a specification), 7 **PARTIALLY COMPLETE** items with a recorded
reason, and 3 **UNVERIFIED** validations (section J, Appendix 1).

**Important remaining risks.**
- No experimental validation; heuristics (primer structure severity,
  guide ranking, hairpin screen) are labeled but remain heuristics.
- Variant normalization is not compared with `bcftools norm`, ORF calls
  not with ORFfinder.
- Remote services (NCBI, literature providers) can change; the live checks
  are snapshots of 2026-11-12/13.
- Dialogs behind files and network (exports, downloads) are rendered but
  only partly driven (ui line coverage 53 %).

---

## B. Scientific improvements

The P0 defects of the previous cycle (1-9) were re-inspected and keep 22
known-answer tests (`tests/test_p0_regressions.py`). Phase 6 corrections:

| Item | Module | Previous behavior | Corrected behavior | Reference / reasoning | Test evidence | Independent validation |
|---|---|---|---|---|---|---|
| D1 | ORF/protein | Initiator Met for GTG/TTG under any code | Only initiation codons of the selected NCBI table | NCBI genetic codes | tier 1 tests | Biopython translation |
| D2 | ORF/protein/RNA | Translation fell back to table 1 | ORF carries its own code | NCBI tables | tier 1 | Biopython |
| D3, D5 | Restriction/cloning | Circular topology ignored | Origin-spanning sites, circular vector/insert | – | tier 1 | Biopython digests, pUC19 |
| D4, P0-5 | Restriction/cloning | Ambiguous sites could cut / define ends | Only CONFIRMED sites cut | IUPAC | tier 1 | – |
| D8-D10, S23, S3-S5 | CRISPR | No cut site; own site counted as off-target; config ignored; Cas12a seed wrong | Modeled cut with source; own site excluded; config honored; PAM-proximal seed | Cas9/Cas12a literature (METHODS) | tier 1/2 | – |
| D12 | Variants | Contig mismatch skipped silently | Flagged | VCF 4.x | tier 1 | – |
| D13/D14 | Input | Peptide in IUPAC letters loaded as DNA | Classifier with evidence levels; user declares ambiguous records | IUPAC | sequence-type tests | Biopython IUPAC sets |
| D14-F6 | Protein | B/Z/J/X/U/O lost or shifted | Kept in place | – | stage 5 tests | ProtParam (31 proteins) |
| D16 | Primer | Shadowed generations | Removed after reachability review, bytecode-identical survivors | – | stage 4 | – |
| D18 | All | Heuristic scores called "confidence" | `heuristic_score`, evidence levels | – | golden (documented change) | – |
| S15-S19 | Variants/cloning/ORF | Bridge broken; digest cut at site start; fixed start set | Fixed | – | tier 2 | – |
| IG weight | Protein | Instability weight IG −7.49 | 1.0 | Biopython DIWV (secondary source) | reference test | Biopython |
| P1-9 | ORF | No origin-spanning or 5'-partial ORFs | Both, with rotation invariance | ORFfinder conventions | `test_orf_partial_circular` | not compared with ORFfinder |
| P2-5 | RNA | G·U pairs ignored | Counted separately | – | `test_rna_wobble` | – |
| P2-7 | Restriction | No Type IIS; forward strand only | BsaI, BsmBI, BbsI, SapI on both strands | REBASE | `test_type_iis` | Biopython 140 cuts; REBASE 611 |
| HGVS | Mutation | `p.R2C`, "Ter" in one-letter code | `p.(R2C)`, `p.(W24*)`, `p.(M1?)`, `p.(*110Qext*?)` | den Dunnen et al. 2016 | `test_hgvs_notation` | no HGVS validator |
| D20 | Primer | Structural penalty of the candidate screen always 0 (read missing fields) | Stems read from the candidate | documented intent of the screen | `test_primer_corrections` | – |
| D21 | Primer | `nearest-neighbor` passed validation, then crashed | One alias set | – | same | – |
| D22 | Mutation | ORF without sequence -> AttributeError | MutationAnalysisError | – | stage 3 test | – |
| D23 | Primer export | Pair risk flags and warnings never exported | Dataclass fields | Master Prompt 6.4 | export test | – |
| D24 | BLAST | Text export and query file written with CR CR LF | Written as built | – | export test | – |
| D25 | BLAST | One network error during polling lost a submitted search and its RID | Bounded retries; RID in the error | found in a live search | `test_blast_polling` | live search succeeded afterwards |

Validation against external references (all reproducible from the suite
unless noted): enzyme definitions 26/26 against REBASE 611 (2026-11-13,
one-off) and against Biopython 1.88; 132 pUC19 digests and 132
pseudo-random digests, 0 differences; end compatibility 484/484; primer Tm
against primer3-py 2.4.0, 2,916 conditions, max |ΔTm| 0.0032 °C; primer
structure screens vs primer3-py ΔG (ρ −0.78/−0.74/−0.51, heuristic);
protein descriptors vs ProtParam (31 proteins) and ExPASy (ubiquitin);
genetic codes vs Biopython; Primer3 manual worked example.

---

## C. Functional completeness

- **Inventory.** `docs/FEATURE_STATUS.md` lists every module, menu entry,
  function group and export with implementation, execution, integration
  and validation status. Every feature of the initial inventory is still
  present; none disappeared.
- **Implemented in Phase 6.** Type IIS digests; circular and 5'-partial
  ORFs; primers around a restriction site with diagnostic digest; variant
  check in primer and CRISPR menus; BLAST RID resume; API entry points for
  every analysis; G·U hairpin count; pipeline parameters for every stored
  result.
- **Corrected.** D1-D25 and the Stage F faults (single-enzyme crash, Enter
  regressions, BLAST hit selector, BLAST query-type question).
- **Tested.** Every menu is rendered and driven (rendered scenarios); the
  protein and RNA views, the export functions and the variant helpers are
  executed with output checks (`test_execution_coverage`).
- **Independently validated.** Restriction definitions and digests,
  primer Tm, protein descriptors, genetic codes (section B).
- **Limitations.** AI Experiment Planner not implemented (labeled);
  off-target search only in supplied sequences; no CFD score; no Golden
  Gate planning; remote BLAST without `core_nt`.

---

## D. Architecture

- **Before.** Twelve analysis modules each held their menus and imported
  `ui.style`; `core.pipeline` imported the interface at call time; the
  pipeline's ORF selection called a dialog.
- **After.**
  ```
  main.py          session; imports dialogs from ui/
  ui/              style.py (renderers), menu.py, report.py,
                   <module>_menu.py (11 dialog modules), variant_check.py
  core/            api.py (every analysis), pipeline.py, models.py,
                   alphabets.py, seqtools.py, orf_context.py, report.py,
                   textformat.py (plain-text layout)
  analyses/        science only; no ui import outside script blocks
  fileio/, data/   readers, enzyme and code tables
  ```
- **Shared components.** `core.textformat` (rules, field lines, columns,
  header text) used by reports and re-exported by `ui.style`; one input
  contract (`_sequence_of`); central alphabets; one restriction matcher;
  one Tm-model alias set; `core.api.binding_features`.
- **Removed duplication.** Superseded primer generations (D16); duplicated
  string layout between reports and screens.
- **API changes.** New: `design_primers_for_restriction_site`,
  `cross_check_variants`, `binding_features`, `blast_resume`. Changed: the
  dialog functions are imported from `ui.<module>_menu`;
  `SynapsePipeline.select_orf(chooser=...)` replaces `interactive=True`;
  `PipelineContext.print_status(render)` takes a renderer.
- **Pipeline changes.** Parameters recorded for every module, including
  the sequence overview, cloning and the restriction-site primer design;
  mapped VCF variants kept per sequence (SHA-256 checked).
- **Data model changes.** `AdvancedPrimerPair` and
  `AdvancedPrimerCandidate` carry their assessment as fields (D23).
- **Modified files.** 42 application files changed since the baseline
  `e4395ed` (18,820 insertions, 16,678 deletions, most of it the P2-2 move).
- **Not done.** Large modules remain large (`analyses/blast.py` about
  4,600 lines, `mutation.py`, `primer.py`); splitting them was not
  required for correctness and was not started.

---

## E. Terminal UI standardization

All eleven module dialogs, the session screens and the pipeline status use
`ui/style.py`: one rule width (70), one meaning of `[0]` (Back/Cancel/
Exit), one prefix per message category, aligned fields, one menu layout
with module-specific groups (Master Prompt 5.3). Per-module history and
the before/after structure are in `docs/STAGE_F_UI_REVIEW.md` and
`docs/STAGE_F_FINAL_REVIEW.md`; the standard and its enforcement rules are
in `docs/TERMINAL_UI_STANDARD.md` (§7 rules run in the suite).

Final-run changes: the dialogs moved to `ui/` with byte-identical
screens; the BLAST hit selector and query-type question follow the
navigation rules; the last hand-padded labels use `style.field`; the
literature limitation is a NOTE.

Remaining inconsistencies (TERMINAL_UI_STANDARD §9): the "RESULT:" lines
of module self-tests run as scripts; plain status lines such as "Loaded 2
FASTA record(s)."; option lists that keep their own prompt label; the
undetermined BLAST query type explained by the raw exception text.

---

## F. Scientific output standardization

Results carry module, data, parameters, warnings, limitations and the
input's SHA-256 (`AnalysisResult`); exports (JSON, Markdown) include
parameters, warnings, provenance and version information (Stage G, 6.4),
and since D23 the primer pair assessment. Units have one notation;
positions shown to users are 1-based and inclusive, API coordinates are
documented per structure; heuristic scores are named as such (D18);
evidence levels for sequence types. Restriction records use one coordinate
convention (`1ae2833`). Remaining: literature dates are kept as each
provider reports them (documented); BLAST E-values are printed with three
significant figures.

---

## G. Language and editorial improvements

The editorial pass of Stage H (2026-11-02) is documented in
`docs/CHANGES.md` (original wording, revised wording, reason, category)
and enforced by `SpellingStandardTests` (American spelling across code,
tests and documentation) and the message-prefix rules. `docs/STYLE_GUIDE.md`
holds terminology, writing, message and code conventions. Text written in
the final run follows it (checked by the same tests). Shared warning and
error conventions: `ERROR:`, `WARNING:`, `NOTE:`, `OK:` through
`ui.style`.

---

## H. Testing

| Item | Result (2026-11-13, HEAD before this report) |
|---|---|
| Total tests | 463 (451 test methods + 12 generated golden-master scenarios) |
| Passed / failed | 463 passed, 0 failed on Python 3.14.6 (1 skipped: the opt-in 10 kb benchmark) |
| Python 3.10.11 | 463 run, 0 failed, 17 skipped (Biopython and primer3-py comparisons, packages not installed for 3.10; one 3.12-tokenizer test; the benchmark) |
| Golden master | 12/12 unchanged in the final run; one documented change in Stage G (D18) |
| Integration tests | `test_integration_workflows` (10.1-10.4, 10.7; 10.2 with restriction-site primers and parameter recording) |
| Rendered interface | `test_rendered_menus` and `test_terminal_ui`: 64 tests, every menu call site rendered |
| Architecture | `LayeringTests`: core imports no ui; every analysis runs with ui unimportable |
| Self-tests | all module self-tests pass inside the suite |
| Reference comparisons | section B; skipped with the reason when a reference package cannot be loaded (Windows application control blocked primer3-py temporarily on 2026-11-13) |
| Benchmark | 10 kb primer design 9.62 s CPU time (machine under background load) |
| Coverage | 71 % combined line/branch (coverage.py 7.10.7; `docs/TEST_COVERAGE.md`) |
| Static analysis | pyright 1.1.414 at Python 3.10: 58 errors; baseline `ab4bd2a` 134; 0 new, 76 removed. Ruff 0.14.4 pyflakes rules: all checks pass after removing 40 unused or duplicated imports; no undefined name was found |
| Environment | Windows 11 Home 10.0.26200; CPython 3.14.6 (`.venv`), CPython 3.10.11 embeddable (`.tools/`); biopython 1.88, primer3-py 2.4.0, numpy 2.5.4, coverage 7.10.7; BLAST+ 2.17.0 (`docs/TOOL_AND_DATA_LEDGER.md`) |

No test was counted as passed without running. Python 3.11-3.13 and
non-Windows platforms were not run.

---

## I. Git history

- 140 commits; 131 in Phase 6 (from `824ca5f`), 124 of them with a
  `Co-Authored-By` line. No history rewrite, no reset, no force-push, no
  remote push.
- By stage: baseline and audit (`e4395ed`-`824ca5f`); Stage B-E corrections
  (2026-10-31); Stage F (2026-11-01 to 11-03); Stages G-I (2026-11-02);
  completion run packages A/B (2026-11-11, `68e59d1`-`15e2ff1`); second
  part (2026-11-12/13, `bf4d096`-the commit adding this report). The full
  list with subjects is Appendix A of
  `docs/PROJECT_SYNAPSE_COMPLETE_HANDOVER.md` and `git log`.
- Since the previous final report (`12949f3`): 30 commits, 68 files,
  16,050 insertions, 9,531 deletions.
- Working tree: clean after the commit of this report; no uncommitted
  changes.

---

## J. Remaining work

| Issue | Module | Why unresolved | Needed | Next action | Kind | Status |
|---|---|---|---|---|---|---|
| `bcftools norm` comparison | variants | No native Windows build | WSL or a Linux host (system-wide: owner approval) | Run the normalization corpus through both | environmental | BLOCKED |
| ORFfinder comparison | ORF | Same | WSL or Linux | Compare ORF calls on reference sequences | environmental | BLOCKED |
| NEBcutter comparison | restriction | Web application without a programmatic interface | Manual use by the owner | Substance covered by REBASE 611 and Biopython | environmental | BLOCKED (substance covered) |
| CFD off-target score | CRISPR | The published matrix is not in the project; its redistribution terms are not established | Matrix file with citation and licence | Implement only with the cited matrix (P2-3) | scientific | BLOCKED |
| Genome-wide specificity | CRISPR, primer | Needs genome-scale data and indexes | Very large downloads: owner approval | P3-2 | environmental | BLOCKED |
| Guruprasad 1990 primary table | protein | The table is in the paper, not machine-readable | The paper | Check the 400 weights against it | scientific | BLOCKED |
| RNAfold | RNA | External binary not installed (system install) | ViennaRNA and owner approval | Optional energy model (P2-5) | technical | BLOCKED |
| Local BLAST on a real database | BLAST | swissprot about 0.5 GB | Owner approval | Run blastp against swissprot | environmental | BLOCKED |
| AI Experiment Planner | main menu [11] | No specification of what it should do; P3-1 deferred | Owner's specification | Labeled "not implemented" | scope | BLOCKED |
| Licence, legacy BLAST result file | project | Owner decisions | Owner | Choose a licence; keep or delete `blast_results/TEST_8_REFERENCE_blast.json` | administrative | BLOCKED |
| Primer generation names (P2-1) | primer | Decided against renaming (no behavior gain) | – | Revisit if the module is split | technical | PARTIALLY COMPLETE |
| `complexity` field name, RNAfold (P2-5) | sequence, RNA | Field kept for compatibility | – | Rename with a deprecation path | technical | PARTIALLY COMPLETE |
| Large modules (7.5) | blast, mutation, primer | Not required for correctness | – | Split by responsibility | technical | PARTIALLY COMPLETE |
| UI §9 remainders | ui | Minor | – | Self-test RESULT lines, status lines | technical | PARTIALLY COMPLETE |
| Dialogs behind files and network | ui | Need files, tools or network | – | Drive export and download dialogs with temporary folders and mocks | technical | PARTIALLY COMPLETE |
| `core_nt` remote database | BLAST | Not in the database list | One live check | Add and verify | technical | PARTIALLY COMPLETE |
| Unrelated long alignments O(n·m) | mutation | Inherent to global alignment | – | Document (done) | scientific | PARTIALLY COMPLETE |
| HGVS validator | mutation, variants | No validator run | Mutalyzer or `hgvs` with its data | Validate a notation corpus | scientific | UNVERIFIED |
| Python 3.11-3.13, Linux, macOS | all | Not run | Interpreters or hosts | CI matrix | environmental | UNVERIFIED |
| Experimental validation | all | Out of scope for software | Laboratory work | – | scientific | UNVERIFIED |

---

## K. Final quality assessment

- **Scientific reliability: good for a screening tool.** Core
  calculations agree with established implementations where compared;
  heuristics are labeled; six defects found in the final run show that
  the test-driven reviews still find real faults, and all six are fixed.
  Not validated experimentally.
- **Functional completeness: high.** Every advertised feature is
  implemented or labeled as not implemented; every menu is rendered and
  driven; workflows 10.1-10.4 and 10.7 are integration-tested, 10.5 and
  10.6 checked live once.
- **Architecture: clean layering.** Science, API, pipeline and interface
  are separated and the separation is tested. Large modules remain.
- **Maintainability: moderate to good.** 463 tests, golden master,
  documented conventions; pyright errors more than halved; long modules
  and historical names (`_v5`-`_v11`) remain.
- **Terminal interface consistency: high.** One renderer and navigation
  model, enforced by tests; minor remainders listed in E.
- **Language quality: consistent.** Style guide enforced by tests.
- **Reproducibility: good.** Parameters, hashes and versions in every
  result; tool and data ledger with checksums; two Python versions.
- **Testing maturity: good.** Unit, regression, golden, integration,
  rendered-interface, architecture and external-reference tests; 71 %
  coverage with the uncovered parts named.

### Independent final audit

The final audit (2026-11-13) re-read the Master Prompt section by section
and checked each requirement against the code, the tests and this
repository's documents (Appendix 1). It re-ran: the full suite on Python
3.14.6 and 3.10.11, the golden master, the self-tests, the rendered and
layering tests, pyright against the baseline, `git diff --check` over every
change since `12949f3`, coverage, the 10 kb benchmark, the REBASE
comparison, and a Ruff pyflakes pass. Findings of the audit itself, all
resolved before this report: two hand-padded labels converted; the pUC19
fixture's trailing blank line kept verbatim (its checksum is tested) and
marked in `.gitattributes`; the REBASE re-check closed open validation 5
of AUDIT_STATUS; the tool evaluation of Master Prompt 7.4 had not been
recorded and was done (Ruff found no undefined name; 40 unused or
duplicated imports, three of them left by the P2-2 move tool, were
removed with screens and suite unchanged).

**Verdict: Phase 6 is PARTIALLY COMPLETE.** Every requirement that could
be completed with the resources and permissions available is COMPLETE and
verified; the items in section J are open for the reasons stated, and none
of them is presented as done.

---

## Appendix 1. Requirement matrix (Master Prompt §0-§18)

| § | Requirement | Status | Evidence |
|---|---|---|---|
| 0 | Coherent, scientifically reliable application; no frontend | COMPLETE | Architecture (D), boundaries (§14) |
| 1.1 | Work only in the copy; inspect before changing | COMPLETE | Appendix 2; initial audit §1 |
| 1.2 | BLAST database policy | COMPLETE | `data/blast_db/` excluded; local tests on small databases; live remote run |
| 1.3 | External research, downloads, traceability | COMPLETE | `docs/TOOL_AND_DATA_LEDGER.md` |
| 1.3.1 | Skills and MCP: availability vs use recorded | COMPLETE | Initial audit; ledger (PubMed tool used once, no skill workflow executed in the final run) |
| 1.4 | Preserve work, history, tests | COMPLETE | No rewrite; golden master; history intact |
| 2.1-2.5 | Full audit and audit report | COMPLETE | `docs/PHASE6_INITIAL_AUDIT.md` |
| 3.1 | P0 re-inspection with reproduction tests | COMPLETE | `test_p0_regressions` |
| 3.2 | ORF analysis | COMPLETE (ORFfinder comparison BLOCKED) | P1-9, D1, D2, D19, S19 |
| 3.3 | RNA: selected ORF vs whole sequence, no silent switch | COMPLETE | pipeline rules, workflow 10.1 test, G·U count |
| 3.4 | Protein translation and analysis | COMPLETE (Guruprasad primary table BLOCKED) | D14-F6, IG weight, ProtParam |
| 3.5 | Primer design, Primer3 comparison | COMPLETE | Tm vs primer3-py; D20, D21, D23; restriction-site design |
| 3.6 | Restriction incl. Type IIS, REBASE check | COMPLETE (NEBcutter BLOCKED, substance covered) | REBASE 611, Biopython, pUC19 |
| 3.7 | Mutation analysis, HGVS distinction | COMPLETE (validator UNVERIFIED) | HGVS-style notation labeled; D22 |
| 3.8 | Variant analysis, bcftools comparison | PARTIALLY COMPLETE (bcftools BLOCKED) | VCF tests, D12, S15/S16 |
| 3.9 | CRISPR incl. P2 items and CFD | PARTIALLY COMPLETE (CFD, genome-wide BLOCKED) | D8-D10, S3-S5, S23, variant check |
| 3.10 | BLAST | COMPLETE (large local database BLOCKED) | coverage P0-1, RID resume, D24, D25, live run |
| 3.11 | Literature | COMPLETE | throttling, retries, live check |
| 3.12-3.13 | Remaining modules; completion requirements | COMPLETE | FEATURE_STATUS |
| 4.1-4.3 | Feature inventory; execution testing; no disguised placeholders | COMPLETE | FEATURE_STATUS; planner labeled |
| 5.1-5.6 | Unified terminal UI and standard | COMPLETE (minor remainders, E) | TERMINAL_UI_STANDARD; rendered tests |
| 6.1-6.4 | Output standardization | COMPLETE | Stage G; D23 |
| 7.1-7.3 | Core architecture, separation, shared functions | COMPLETE | P2-2; LayeringTests |
| 7.4 | Engineering standards, static analysis | COMPLETE | pyright 134 -> 58; Ruff pyflakes pass clean; Black, mypy, pytest evaluated (`docs/TOOL_AND_DATA_LEDGER.md`) |
| 7.5 | Maintainability (no giant modules) | PARTIALLY COMPLETE | large modules remain |
| 8 | Language and editorial standardization | COMPLETE | Stage H; spelling tests |
| 9 | STYLE_GUIDE | COMPLETE | `docs/STYLE_GUIDE.md` |
| 10.1-10.4, 10.7 | Workflows | COMPLETE | integration tests |
| 10.5, 10.6 | BLAST and literature workflows | COMPLETE | parsing tests, live checks |
| 10.8 | Integration requirements | COMPLETE | stale-result invalidation, parameter recording, SHA-256 checks |
| 11.1-11.7 | Test levels, golden protection | COMPLETE | section H |
| 11.8 | Reproducibility across Python versions | PARTIALLY COMPLETE | 3.10 and 3.14 run; 3.11-3.13 UNVERIFIED |
| 11.9 | Full suite after each stage | COMPLETE | each commit message |
| 12 | Feature matrix kept current | COMPLETE | FEATURE_STATUS (2026-11-13) |
| 13 | Staged, controlled implementation | COMPLETE | PHASE6_PROGRESS rows 1-60 |
| 14.1-14.6 | Boundaries | COMPLETE | no frontend, no fabrication, approvals respected |
| 15 | Documentation | COMPLETE | README, CHANGES, METHODS, AUDIT_STATUS, FEATURE_STATUS, standards |
| 16 | Final criteria | PARTIALLY COMPLETE | section J |
| 17 | Final report A-K | COMPLETE | this document |
| 18 | Final directive | PARTIALLY COMPLETE | open items stated, not hidden |

## Appendix 2. Project integrity

- Work only in `C:\Users\bayra\OneDrive\Desktop\Project Snyapse - Kopie`.
  The original folder `Project Synapse` was only read (the Master Prompt
  PDF); nothing was written there.
- ReadOnly attribute: 15 of 15 protected files set at every commit of the
  final run (`main.py`, ten analysis modules, `core/pipeline.py`, two data
  modules, `ui/menu.py`); cleared only for the files a package edited.
- Git: no history rewrite, reset, force-push or remote push; every package
  committed separately; working tree clean.
- Approvals: installations and downloads followed the policy the owner
  adopted on 2026-11-12 (project-local packages and interpreters, small
  public datasets); no project file, private sequence or personal data was
  sent to an external service; a Windows application-control block of a
  compiled test dependency was not bypassed.
- No scientific data, citation or validation result was invented; every
  comparison named here was run, and the ones not run are listed in J.
