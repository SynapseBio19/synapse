# Project history and decisions

A concise, traceable record of how Project Synapse reached its current
state, which decisions bind further work, and what is still open. It
complements `docs/PROJECT_SYNAPSE_COMPLETE_HANDOVER.md` (the long dossier)
and `docs/PHASE6_FINAL_REPORT.md` (the current verdict).

Written 2026-11-13 at the start of the continuation run that follows the
owner's "Master Execution Prompt" (full project continuation, historical
chat analysis, Phase 6 completion and final audit).

---

## 1. Sources and their limits

| Source | What it covers | Analyzed |
|---|---|---|
| Claude Code session "Project Synapse Phase 6" (`local_3b57cea9`, 2026-10-30 to 2026-11-11, about 1,600 messages, several context summaries) | Phase 6 audit, Stages B-I, Stage F re-examination, completion run packages A and B, handover dossier | Yes: all 47 owner messages extracted from the exported transcript and read; assistant replies read at the decision points (approvals, Stage F review, the literature-search fault) |
| Claude Code session "Project Synapse Handover Setup" (`local_c073d91f`, 2026-11-11 to 2026-11-13) | Approval policy, completion run T1-T12, D20-D25, final report A-K, comparison with the Master Prompt PDF | Yes: all 11 owner messages and the closing answers |
| Claude Code sessions "Reload skills", "Available plugins, skills, hooks, and tools" (2026-10-30) | Tool setup only | Titles and scope checked; no project decisions |
| Conversations in the claude.ai Project, and the earlier "Claude Cowork" development cycle | Original vision, the pre-Phase-6 audit, P0 fixes, Stages 0-5 | **Not accessible from this environment.** Reconstructed only from Git (`e4395ed`-`0dad3bb`, 2026-10-29), `docs/AUDIT_STATUS.md`, `docs/CHANGES.md` and Master Prompt §0 |
| Git history (142 commits, 2026-10-29 to 2026-11-13) | What was actually changed | Yes |
| Phase 6 Master Prompt, `Project_Synapse_Phase_6_Master_Prompt_Skills_Integrated.txt` (2,344 lines, SHA-256 `c4a38c83…07346`, in the owner's Downloads folder) | Governing technical specification | Yes, in full. A shorter "Revised" text (2,256 lines) differs only by the skills/MCP section §1.3.1 and five cross-references to it; the PDF in the original folder was compared with the text version in an earlier session |

Rule applied throughout (Master Execution Prompt §1.3): conversations give
intent and history, the repository gives the implementation, tests give
what is verified, and external sources decide scientific claims.

---

## 2. Original vision

- An intelligent bioinformatics platform that supports biologists,
  biotechnology students and researchers in sequence analysis: one
  coherent application in which modules pass selections and results to
  each other, not a set of scripts (Master Prompt §0).
- Its value lies in integration, scientific transparency, reliable
  calculations and an understandable interface. It does not replace
  Benchling, Geneious, AlphaFold, a LIMS, or established tools such as
  Primer3, BLAST, ORFfinder, bcftools, NEBcutter or RNAfold.
- Heuristics must never be presented as validated probabilities.
- Origin: `LEGACY/first_dna_analyzer.py`, a single nucleotide counter.
- A graphical or web front end is a later phase (P3-3). `core.api` is the
  prepared foundation; front-end work is out of scope (§14.1).

## 3. Major development stages

| Period | Stage | Key commits |
|---|---|---|
| 2026-10-29 | Previous cycle (Claude Cowork): baseline, test suite and golden master, nine P0 fixes, shared data model, P1 improvements, core API, METHODS, AUDIT_STATUS; 109 tests | `e4395ed`-`0dad3bb` |
| 2026-10-30 | Phase 6 start: pyright configuration, skills ignored in Git, initial audit D1-D19 / S1-S14 and feature matrix | `5242fcb`, `5d8c73a`, `824ca5f` |
| 2026-10-31 | Stages B-E: D1-D19 corrections one commit each; S-findings; superseded primer generations removed; RID resume; API entry points | `dc780f6`, `ab4bd2a` (pyright baseline 134), `dfe773d`, `ad843c5` |
| 2026-11-01 to 11-03 | Stage F (terminal UI), re-examined after the owner judged it incomplete; 16 interface files converted; core no longer imports ui; every menu rendered and driven | `5db21ab`-`12949f3` |
| 2026-11-02 | Stages G-I: heuristic scores renamed (D18), evidence levels, editorial pass, STYLE_GUIDE, primer Tm vs primer3-py, IG instability weight, first final report | `835744e`, `2851e11`, `1f0c0f5` |
| 2026-11-11 | Completion run packages A and B: Type IIS, HGVS-style notation, circular and 5'-partial ORFs, G·U count, alignment memory and speed, cloning planner cache | `68e59d1`-`15e2ff1` |
| 2026-11-12 to 11-13 | Completion run T1-T12, D20-D25: primer 10 kb in about 10 s, restriction-site primers (10.2), P2-2 dialogs moved to `ui/`, coverage, ledger, Python 3.10 run, final report A-K | `aba3fb0`-`ea9612b` |

## 4. Binding decisions (do not reverse without the owner)

### Owner decisions

| Date | Decision | Source |
|---|---|---|
| 2026-10-30 | Work only in `Project Snyapse - Kopie` (folder spelling kept); never modify the original `Project Synapse` folder | Master Prompt §1.1, every owner prompt |
| 2026-10-30 | Biopython and primer3-py stay strictly optional reference/test dependencies, never runtime dependencies (`dependencies = []`) | owner prompt 2026-10-30 11:14 |
| 2026-10-30 | Keep all eleven local K-Dense skills; ignore exactly `/.claude/skills/` in Git | owner approval 2026-10-30 22:50 |
| 2026-10-30 | ReadOnly attribute of the 15 protected files: clear only on the files a task needs, never recursively; restore after the commit | owner approval 2026-10-30 23:06 |
| 2026-10-30 | Audit documents must not claim states that were not verified (the "clean working tree" correction) | owner approval 2026-10-30 22:50 |
| 2026-11-02 | Stage F is judged by actual visual and structural consistency, not by shared separators and prefixes; module-specific menu structures may differ | owner review request 2026-11-02 18:15 |
| 2026-11-03 | A local import is not architectural decoupling: core must import and run without `ui` (tested) | owner release of recommendations 4 and 5 |
| 2026-11-12 | Approval policy = Master Prompt §1.3/§14.5: project-local installs, small and moderate downloads and public sources are allowed; ask first for very large downloads, system-wide changes, costs, uploads of private data, irreversible actions | owner message 2026-11-11 23:34 ("diese Regeln übernehmen") |
| 2026-11-13 | **AI Experiment Planner stays COMING SOON.** It is not implemented, researched or connected in this execution; the menu entry and placeholder stay | Master Execution Prompt §3 (supersedes the earlier "label as not implemented / needs a specification" handling, see §6) |
| 2026-11-13 | Communicate with the owner in German; all code, comments, output and technical documentation in English | Master Execution Prompt §2.1 |

### Engineering decisions recorded in the repository

- Science lives in `analyses/` and `core/`; dialogs live in `ui/`
  (`50dea66`). `LayeringTests` enforce that core and analyses run with
  `ui` unimportable.
- One terminal renderer, `ui/style.py`; no second UI framework.
- User-facing positions are 1-based inclusive; API coordinates are
  documented per structure; restriction site and cut records share one
  convention (`1ae2833`).
- Pipeline rule: a selected ORF feeds RNA and protein; without a
  selection RNA analyzes the whole selected sequence; the ORF's genetic
  code travels with it (D1, D2).
- Heuristic scores are called `heuristic_score`, never "confidence" (D18).
- Golden-master snapshots change only with a documented scientific
  reason (one change so far, D18).
- Primer generation names `_v5`-`_v11` are kept (P2-1): renaming changes
  no behavior and was deferred until the module is split.
- No history rewrite, reset or force-push; one commit per package with
  the evidence in the message.

## 5. Important corrections that came from the owner

- **Stage F was overstated** (2026-11-02). The owner saw that the
  interface was not unified; the re-examination found that the shared
  renderers existed but no module used them. The conversion of 16 files
  and the rendered-menu tests followed.
- **The literature search returned nothing** (2026-11-03). The Stage F
  conversion had turned the hint "[none]" into a value; six prompts
  (literature dates, primer conditions, local BLAST options, cloning
  strategy count) broke on Enter. Fixed in `1516339` with tests that
  drive the dialogs with empty input and a source rule against the fault.
- **Final reports were too optimistic** twice (2026-11-02, 2026-11-13).
  The owner asked for independent final controls; the reports were
  rewritten with COMPLETE / PARTIALLY COMPLETE / UNVERIFIED / BLOCKED.
- **Context budget.** The owner stopped the run after package B
  (2026-11-11) to continue in a fresh session; handover documents and the
  continuation record exist for that reason.

## 6. Contradictions found and how they are resolved

| Topic | Statements | Resolution |
|---|---|---|
| AI Experiment Planner label | Initial audit D15 (2026-10-30): "coming soon" advertises a feature that does not exist; relabeled "not implemented" (`f3a09c1`). Final report 2026-11-13: BLOCKED, "needs the owner's specification". Master Execution Prompt (2026-11-13): status must be COMING SOON, keep the entry, do not implement | The later explicit owner instruction wins. The entry reads "Coming Soon" and its message says that the feature is planned but not yet available and that nothing is run or stored, which keeps D15's honesty requirement. It is counted neither as a completed module nor as a blocked requirement but as **intentionally deferred** |
| Approval for installs and requests | Handover 2026-11-11: install nothing, no external requests. Master Prompt §1.3 and owner 2026-11-12: routine project-local installs and public sources allowed | Owner decision of 2026-11-12 applies (§4) |
| Historical figures | Earlier records: 109, 228, 267, 283, 314, 326, 421 tests; pyright 134 -> 132 -> 112 -> 111 -> 58 | Each is a snapshot of its date. Current values are measured, not copied (§7) |
| "Phase 6 complete" | Final report of 2026-11-02 declared completion; the independent review the owner requested found it overstated | Superseded; the report of 2026-11-13 says PARTIALLY COMPLETE |
| Primer 10 kb timing | 135.6 s -> 25.6 s (2026-11-11), -> about 10 s / 9.6 s CPU (2026-11-12/13) | Later measurement is current; it is re-measured in this run |

## 7. State at the start of this run (measured 2026-11-13)

- `main` at `ea9612b`, working tree clean, 142 commits.
- 463 tests, all passing on CPython 3.14.6 (1 skipped: the opt-in 10 kb
  benchmark); golden master 12/12 unchanged.
- pyright 1.1.414 at Python 3.10: 58 errors (baseline `ab4bd2a`: 134).
- ReadOnly set on all 15 protected files.
- `.venv`: biopython 1.88, primer3-py 2.4.0, numpy 2.5.4, coverage
  7.10.7, ruff 0.14.4, pyright 1.1.414. Project-local CPython 3.10.11 in
  `.tools/`. WSL is not installed; the system has only Python 3.14.

## 8. Outstanding work as recorded on 2026-11-13

From section J of the final report, re-assessed for this run:

- Validations previously blocked: bcftools normalization, ORFfinder,
  NEBcutter, CFD matrix, genome-scale specificity, the Guruprasad 1990
  primary table, RNAfold, local BLAST on a real database, HGVS validator,
  Python 3.11-3.13.
- Partial items: `core_nt`, file and network dialogs only partly driven,
  UI §9 remainders, the `complexity` field name, large modules, legacy
  primer names, pyright diagnostics.
- Owner decisions still pending: the project licence
  (`pyproject.toml`: "Proprietary - all rights reserved (decision
  pending)") and whether to keep `blast_results/TEST_8_REFERENCE_blast.json`.

What this run did with each item is recorded in `docs/PHASE6_FINAL_REPORT.md`
and `docs/PHASE6_CONTINUATION.md`.

## 9. Context that is not visible in the code

- The owner reviews claims independently and expects evidence: tests that
  fail on the previous code, byte-identical output for refactors, and
  explicit BLOCKED items rather than silent omissions.
- During Phase 6 a CRISPR-related investigation (S3, S4, S5, S23) was
  interrupted by a model safeguard on 2026-10-31 and completed under a
  separately scoped assignment; the findings are in the initial audit.
- `pyright.exe` is blocked by Windows application control; run pyright as
  `python -m pyright`. The same control once blocked primer3-py's compiled
  module (2026-11-13); reference tests skip with that reason instead of
  bypassing the control.
- `Synapse_Website/` exists from an earlier period and is out of scope.
