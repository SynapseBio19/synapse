# Phase 6 closing review

Review of the Phase 6 end state at `1f0c0f5`, requested by the owner on
2026-11-02 before accepting the project state. No source file was changed
for this review, no external request was made, and the original project
was not touched. Every check below was run against the repository; the
commands and their outputs are summarized with each result.

## 0. Scope, method and a limitation

**Questions** (from the owner): (1) Master Prompt requirements met or
documented as open; (2) the last 18 commits against the final report;
(3) the instability-index correction; (4) the golden-master change;
(5) whether the 314 tests cover what is claimed; (6) unintended damage to
functions, interfaces or workflows; (7) completeness and consistency of
the open items.

**Method.** Checks were mechanical wherever possible, so that their result
does not depend on the reviewer's memory of the work:

- **Cross-run of test suites.** Two scratch trees were built with
  `git archive` (outside the repository): the `aa5703b` code with the
  current tests, and the current code with the `aa5703b` tests. Every
  failure in either direction was classified.
- **Interface diff.** The top-level functions, classes, constants, class
  fields and argument lists of every changed source file were compared by
  AST between `aa5703b` and `HEAD`.
- **Golden master.** Both changed snapshots were loaded at `aa5703b` and
  `HEAD`; the key set and the values with the key mapped back were
  compared.
- **Commit audit.** Files touched per commit, word-level diff counts and
  attribution trailers were compared with the claims in the documents.
- **Document consistency.** Open-item lists, counts and status statements
  were compared across `PHASE6_FINAL_REPORT.md`, `PHASE6_HANDOVER.md`,
  `PHASE6_PROGRESS.md`, `AUDIT_STATUS.md` and `FEATURE_STATUS.md`.

**Limitation.** This review was carried out by the same model that did
the Stage G-I work in the same session. The mechanical checks above
reduce, but do not remove, that conflict of interest: judgment calls
(what counts as an adequate editorial review, whether a requirement is
"met") are still the author's. A review by the owner or by a separate
reviewer without this session's context remains advisable.

---

## 1. Results by question

### Q1. Master Prompt requirements

Assessed against the completion criteria of section 16 and the
traceability rules of sections 1.3, 1.3.1 and 11.8.

| Area (section 16) | Verdict | Evidence / gap |
|---|---|---|
| Scientific quality: modules audited, defects corrected, assumptions and validation documented, external comparisons distinguished | **Met, with gaps** | D1-D19 corrected; reference comparisons recorded with their level. Gaps: S9 (alignment memory) is quantified, not corrected; four bibliographic entries added in Stage H/I were written without verification (section 2). |
| Functional quality: every feature accounted for, executed, workflows work, results stored and transferred | **Met** | Feature matrix refreshed; 7 workflow tests; the old suite still passes on the new code except for intended renames (Q6). |
| Architecture: shared logic, reduced duplication, core independent of rendering, API/pipeline consistent | **Partly met, documented** | `core/api.py` complete, `ui/style.py` shared; menus still inside the analysis modules (P2-2); `core/pipeline.py` still prints. |
| Python engineering standards (7.4) | **Partly met** | pyright reports **112 errors** in standard mode (29 in `primer.py`, 11 in `crispr.py`, 10 in `mutation.py`, about 35 in tests; mostly `reportArgumentType` and `reportOptionalMemberAccess`). Down from 134, but not clean, and the reports present the number only as "no new diagnostics". |
| Terminal interface | **Met, deviations documented** | 12 enforcement tests; three documented deviations. |
| Scientific communication, style guide reflects implementation | **Mostly met** | Inventory-based review; enforcement tests. The review of long scientific explanations was by inventory and sampling, not line by line (stated in report K). |
| Testing: suite, integration, regressions, reference comparisons, reproducible environments | **Partly met** | 314/314 pass, 0 skipped. Not met: reproducibility across supported Python versions (only 3.14 available); no line-coverage measurement (no coverage tool installed). |
| Documentation accurate and current | **Not fully met** | Several stale or inaccurate statements (section 2 and Q7). |
| 1.3 / 1.3.1 / 11.8 traceability of tools and MCP use | **Partly met** | Versions of Biopython 1.88, primer3-py 2.4.0, BLAST+ 2.17.0, pyright 1.1.414 and Expat are recorded. There is no consolidated record with source URL, access date, license and checksum per tool or dataset (pUC19 L09137.2, BLAST+ installer), and earlier PubMed MCP use is recorded only as "checked in PubMed (PMIDs ...)" without the eight points of 1.3.1. |

### Q2. The last 18 commits against the final report

**Confirmed.** All 18 commits since `aa5703b` exist as listed in report
section I and handover section 12, each is one package, each carries the
Opus 5.5 trailer, and the files each touched match its described scope.
No commit changed a source file outside its stated purpose.

**Inaccuracies found in the documents** (none affects code):

| Location | Statement | Actual |
|---|---|---|
| `CHANGES.md` (Stage H), `PHASE6_PROGRESS.md` row 35, `PHASE6_FINAL_REPORT.md` G | "92 words in 10 documents and 6 code text sites" | Word-level diff of `1a8f2bc`: 86 words in 8 documents; 9 sites in 4 source files (`crispr.py` 6, `api.py` 1, `loader.py` 1, `orf_context.py` 1) plus 6 in tests. |
| `PHASE6_FINAL_REPORT.md` A | "Stage I found and corrected three further defects (items 25-27)" | Item 27 (primer coordinates) was found and fixed in Stage G (`ebf4d52`); Stage I found two. |
| `PHASE6_FINAL_REPORT.md` A | "109 at the start of the previous cycle" | `AUDIT_STATUS.md` records 109 at the **end** of the 2026-10 cycle. |
| `PHASE6_FINAL_REPORT.md` A | "every S item is assessed: confirmed ones are corrected" | S9 is confirmed and only quantified; S14 is partly resolved (menus remain, P2-2); S5 is partly confirmed and addressed. |
| `AUDIT_STATUS.md` P1-1 | "312 at the end of Phase 6" | 314 (two tests were added afterwards in `db74504`). |
| `PHASE6_FINAL_REPORT.md` header | "Work by Claude Opus 5 and Claude Opus 5.5, as recorded per commit" | 5 of the 63 Phase 6 commits (`824ca5f`, `d4544ad`, `5c098e3`, `47146b0`, `7334296`) carry no model attribution. |
| Commit message `2851e11` | P1-3 criterion "had never been run" | Already corrected in `CHANGES.md`; history unchanged by design. |

### Q3. The instability-index correction

**Implementation: confirmed.**

- The change is one value in `INSTABILITY_ROWS["I"]`, column G (column
  order `ACEDGFIHKMLNQPSRTWVY`): -7.49 -> 1.0. The value -7.49 had been
  present since the baseline commit `e4395ed`.
- Table orientation matches Biopython's documented convention
  (`DIWV['D']['G']` is the dipeptide DG): the other 398 of the 400 weights
  equal Biopython's, and FY differs only slightly (33.6 vs 33.601).
- The formula is 10/L times the sum of dipeptide weights with L = number
  of scored residues, as in Guruprasad et al.
- On the old code, 20 of 31 protein subtests, the table test and the
  known-answer test fail; on the new code all pass.
- Effect size: one IG dipeptide changes the index by 8.49 x 10/L (1.70 at
  L = 50, 0.28 at L = 300), enough to move a protein across the
  conventional stability threshold of 40.
- The claim that the audit's proteins contain no IG holds for ubiquitin
  (checked); for lysozyme it is an inference from the audit's matching
  values, because the sequence is not stored in the repository.

**Scientific evidence: plausible, not sufficient.**

- The only reference is Biopython's `ProtParamData.DIWV`. Its source
  comment cites Guruprasad K., Reddy B.V.B., Pandit M.W., Protein
  Engineering 4:155-161 (1990), which corroborates the citation added to
  METHODS, but not the value.
- The documents state that Biopython's table "reproduces ExPASy
  ProtParam". That statement was not verified in this phase: Biopython
  cites the paper, not ExPASy.
- The "transcription error" explanation (the value -7.49 occurs at IK and
  IV in the same row) is a plausible hypothesis, not evidence.
- So the correction rests on one secondary source. It should be
  confirmed against the published table or against ExPASy ProtParam on a
  protein containing IG before it is described as established. The
  documents already list this as open, but they word the ExPASy link more
  strongly than the evidence allows.

### Q4. The golden-master change

**Confirmed: only the intended rename.**

- In the last 18 commits only `tests/golden/orf.json` and
  `orf_legacy_baseline.json` changed, both in `317097d`.
- Key sets differ only by `confidence` (old) vs `heuristic_score` (new),
  6 and 23 occurrences. With the key mapped back, both files are equal to
  their `aa5703b` versions. SHA-256 prefixes before (`e5fe5dd9`,
  `67fef456`) match those recorded in `CHANGES.md`.
- The insertion/deletion counts (458/458) come from the sorted-keys
  writer.
- In `analyses/orf.py`, the diff of `317097d` reduces, after mapping the
  new names back, to the section heading, two comments, the alias line,
  one docstring word and one interpretation sentence. No logic changed.
- The view key in `tests/golden_scenarios.py` was renamed in the same
  commit and still reads the same score.

### Q5. Do the 314 tests cover what is claimed?

**Confirmed.**

- 314 run, 314 pass, 0 skipped, 0 failed at `1f0c0f5`. BLAST+,
  Biopython and primer3-py are present, so no test was skipped.
- The current tests on the `aa5703b` code: 63 failures, all in tests
  added for a Stage G-I change (D18 vocabulary, notation, primer
  coordinates, exported warnings, restriction parameters, IG weight,
  Python 3.10 f-string, spelling, style guide, golden key). This confirms
  that each regression test detects the defect it was written for.
- The characterization tests pass on both code states, as they should:
  the primer3-py Tm comparison, the restriction and end-compatibility
  references, and five of the seven workflow tests (the two that fail
  on the old code check the restriction parameters added in `81ccc55`).
- Every test named in the style guide's enforcement table exists
  (itself tested).

**Not covered or not measurable.**

- **Line coverage is unknown.** No coverage tool is installed, and
  installing one was outside this review's limits.
- **Edge cases of 11.4.** Covered: empty input, multiallelic VCF
  (`A,T` split), unknown file format, circular and origin-spanning
  sequences, minus strand, partial ORFs, missing contigs. Not found as
  explicit tests: a missing or unreadable input file, and an invalid menu
  selection fed to a running menu (the interface tests check labels, not
  this runtime path).
- The 10.5 BLAST and 10.6 literature workflows are covered offline or
  with stubs only.

### Q6. Unintended damage to functions, interfaces or workflows

**No unintended damage found.**

- The `aa5703b` suite on the current code: 283 run, 23 fail, every one
  traceable to an intended change. 22 are D18 renames
  (`SequenceClassification.confidence`, `SequenceRecord.type_confidence`,
  the `confidence` and `sequence_type_confidence` keys); 1 is the exit
  message ("Thank you for using Synapse!" -> "...Synapse."). Every other
  old test passes, including the end-to-end menu runs, the golden master
  apart from the renamed key, and all module self-tests.
- AST interface diff `aa5703b..HEAD`:
  - Removed without an alias: `core.alphabets.CONFIDENCE_LEVELS`,
    `SequenceClassification.confidence`, `SequenceRecord.type_confidence`.
    All three were introduced in Phase 6 (`d057ac1`, `776de4a`), so no
    pre-Phase-6 caller can depend on them. They are documented in
    `CHANGES.md`.
  - `calculate_orf_confidence` remains callable as an alias.
  - Added: `calculate_orf_heuristic_score`, `EVIDENCE_LEVELS`,
    `SequenceClassification.evidence`, `SequenceRecord.type_evidence`,
    `embedded_warnings`, `AnalysisResult.__post_init__`, `_evalue_text`,
    `_enzyme_database_version`.
  - No argument list of a public function changed.
- **Intended behavior changes a user or caller will notice** (all
  documented):
  - Exported JSON keys renamed (D18).
  - The primer menu now reads 1-based inclusive positions. A user who
    types the old 0-based values gets a target shifted by one base.
  - Report summaries now contain payload warnings.
  - BLAST text reports show E-values with three significant figures.

### Q7. Open items: complete and consistent?

**Inconsistencies.** Three documents still contain sections that
contradict the end state:

- `PHASE6_HANDOVER.md` opens with "Phase 6 is **not complete**. Stages
  F-H ... and Stage I ... have not been started", and its section 1
  table gives HEAD `f803577` and 283 tests. The later addenda correct
  this, but a reader starting at the top is misled.
- `PHASE6_PROGRESS.md` "Paused (2026-10-31)" still says "Stages F-I have
  not started", and "Open after the current package" still lists D17 and
  D18 as open.
- `PHASE6_PROGRESS.md` classification row P2-5 still says the complexity
  renaming "belongs with Stage G terminology"; it was resolved by
  documentation in `1269335`.

**Open items listed elsewhere but missing from report section J:**

- G·U pairing and the optional RNAfold integration (P2-5).
- Merging the surviving primer generations (P2-1).
- The NEBcutter comparison (P1-4; `FEATURE_STATUS.md` calls it superseded,
  `AUDIT_STATUS.md` still lists it as open).
- A local BLAST run on a real database such as swissprot.
- An HGVS conformance check against a validator (J lists only the
  parentheses convention).
- `core/pipeline.py` printing from `core/` (a follow-up in
  `TERMINAL_UI_STANDARD.md`).
- The "complexity" field names kept for compatibility.
- The test-hygiene `ResourceWarning` (mentioned in H, not in J).
- The traceability gaps of 1.3.1 / 11.8.
- The 112 remaining pyright errors.

---

## 2. Summary

### Confirmed results

1. The 18 commits after `aa5703b` match their documented scope. Each is
   one package, attributed to Opus 5.5, and touches no source file
   outside its purpose.
2. The golden-master change is exactly the `confidence` ->
   `heuristic_score` key rename (6 + 23 keys). All values are identical,
   and the recorded hashes match.
3. The IG correction is implemented as described: one weight, correct
   orientation, correct formula. Its tests fail on the old code and pass
   on the new.
4. All 314 tests pass, none skipped. Every regression test added in
   Stages G-I fails on the code before its fix.
5. No unintended breakage. Every old-test failure on the new code is an
   intended, documented rename or wording change. No public function
   changed its arguments, and the removed names are Phase 6's own.
6. D18 changed names only; ORF logic is unchanged.

### Claims not sufficiently substantiated

1. The IG weight's correctness rests on one secondary source (Biopython).
   "Reproduces ExPASy ProtParam" is unverified, and the paper's table was
   not consulted.
2. "The CRISPR module was unimportable on Python 3.10/3.11." This follows
   from the pre-PEP-701 grammar, but it was not observed on a 3.10 or
   3.11 interpreter.
3. Four references added to `METHODS.md` in this session (Gardiner-Garden
   and Frommer 1987; Guruprasad et al. 1990; Cock et al. 2009;
   Untergasser et al. 2012) were written without verification. Only the
   Guruprasad journal, volume and pages are corroborated, by Biopython's
   source comment.
4. The editorial review of long scientific explanations was by inventory
   and sampling. The report states this, but "existing text has been
   reviewed" (section 16) is met only in that sense.
5. The counts and phrasings listed in Q2 (92 words, three Stage I defects,
   109 "at the start", "confirmed ones are corrected", 312 tests, "as
   recorded per commit").
6. "Lysozyme contains no IG" is inferred from matching audit values; the
   sequence is not in the repository.

### Open risks

1. **Python version risk.** No run on the minimum supported Python 3.10.
   A static check found one blocker; others of a kind the check does not
   cover (runtime APIs, behavior differences) cannot be excluded.
2. **Type safety.** 112 pyright errors remain, concentrated in `primer.py`,
   `crispr.py` and `mutation.py`. Some `reportOptionalMemberAccess`
   findings may be reachable `None` dereferences.
3. **Scientific references not confirmed at the source**: the instability
   table, the REBASE enzyme data (P0-4), and the four unverified citations.
4. **External services untested live**: literature (S11), remote BLAST
   and RID resume.
5. **Interface change for users**: the primer menu's switch to 1-based
   input is documented but silent at run time beyond the prompt text.
6. **Unknown line coverage**: test breadth is good, but untested code paths
   cannot be quantified.
7. **Documentation drift**: the stale sections in Q7 can mislead the next
   developer about the state of the project.

### Concrete next steps

These are proposals. None was carried out, and those marked "approval"
need the owner's decision.

1. **Correct the documents** (documentation only, no code): the Q2
   inaccuracies, the stale sections in the handover and progress records,
   the P2-5 row, and the missing items in report J. Add a status line at
   the top of `PHASE6_HANDOVER.md`.
2. **Mark the four unverified references** as unverified in METHODS, or
   verify them (approval: PubMed).
3. **Confirm the IG weight** against the Guruprasad table or ExPASy
   ProtParam on a protein containing IG (approval: external access). Until
   then, word the ExPASy link as "Biopython's table, which cites
   Guruprasad et al. (1990)".
4. **Run the suite on Python 3.10** (approval: install a project-local
   3.10 interpreter).
5. **Measure line coverage** (approval: install `coverage` into `.venv`;
   small, test-only).
6. **Add a 1.3.1 / 11.8 tool and data ledger**: name, source URL,
   version, access date, license, checksum, how it was used and what it
   showed. Cover Biopython, primer3-py, BLAST+, pyright, pUC19 L09137.2
   and the PubMed MCP lookups.
7. **Triage the 112 pyright errors** in the three analysis modules for
   real `None` paths; add tests where a path is reachable.
8. **Add the two missing edge-case tests**: a missing input file and an
   invalid menu selection at run time.
9. **Make the primer input change visible**, for example with a one-time
   note in the primer menu ("positions are now 1-based"), if users may
   rely on the old convention.
10. **Ask an independent reviewer** to repeat Q1 and Q3 without this
    session's context.
