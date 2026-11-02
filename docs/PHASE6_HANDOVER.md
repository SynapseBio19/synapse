# Phase 6 handover – planned pause (2026-10-31)

This is the state at a planned pause requested by the owner. Phase 6 is
**not complete**. Stages F–H (terminal UI standard, output
standardization, style guide/editorial pass) and Stage I (final
verification and report) have not been started. The running record is
`docs/PHASE6_PROGRESS.md`; the findings are in
`docs/PHASE6_INITIAL_AUDIT.md`.

## 1. Repository state

| Item | State |
|---|---|
| Branch | `main` |
| HEAD | `f803577` (see section 10; `642f31b` was the state at the original pause) |
| Working tree | clean |
| Tests | 283/283 pass, 0 skipped (`python -m unittest discover -s tests -t .`, ~26 s) |
| Golden master | unchanged since the audit; no snapshot was regenerated |
| pyright 1.1.414 (`python -m pyright --pythonversion 3.10`) | 112 errors; baseline at `ab4bd2a` was 134; no new diagnostic in any package. 22 of the removals came with the dead primer code in D16 part 2. `pyright.exe` is blocked by a Windows application control policy (WinError 4551); use `python -m pyright`. |
| ReadOnly attributes | restored on all 15 files (see the handover notes in `PHASE6_PROGRESS.md`) |
| Optional/test-only tools present | Biopython 1.88, primer3-py 2.4.0, BLAST+ 2.17.0; bcftools absent |

## 2. A note on the pause instructions

Two pause instructions arrived while work was running. The second ("stop
after D8; do not start D10, D11, D15, D16, D19 or D17/D18") reached the
session only after D8 *and* the following packages had already been
completed: D10/D11, D19, D15, S15/S16, S17 and S19. Each of these is a
separate, tested, documented commit (listed below). Nothing was reverted,
because reverting validated work would be a destructive change; the
owner can still revert any single commit with `git revert <hash>`. No
work was started after the instruction was seen: D16, D17, D18, the
remaining S-findings and Stages F–I are untouched.

## 3. Commits of this execution (oldest first, after `ab4bd2a`)

| Commit | Package | Validation |
|---|---|---|
| `02d2899` | Progress record; D14 sub-findings F1–F12, S20–S22 recorded | docs |
| `b8f5ea8` | D14-F2 one BLAST query-type vocabulary (protein BLAST was unusable) | internal tests; external tool execution (makeblastdb + blastp) |
| `d057ac1` | Central `core/alphabets.py` classifier (CONFIRMED / PROBABLE / AMBIGUOUS / INVALID, no frequency threshold) | internal tests; Biopython `IUPACData`; PubMed-checked sources |
| `776de4a` | D14-F3, F11 molecule type from GenBank/FASTA; type in `SequenceRecord`, provenance, report | internal tests; NCBI `[moltype=]` documentation |
| `7c716a2` | D14-F1, S20 BLAST uses the classifier; J/U/O accepted | internal tests; blastp run; BLAST+ U/O/J acceptance checked |
| `162ca16` | D13 + D14-F4/F5/F7 one input contract; menu type gate with user declaration; API type checks; protein-record route | internal tests; end-to-end `main.py` |
| `343ffb7` | D14-F6 protein codes keep position; Sec/Pyl masses; B/Z/J ranges | internal tests; Biopython masses; one Stage 5 assertion updated (it encoded the deletion) |
| `d40d328` | D14-F8/F9/F10 ORF search, transcription, GC counting | internal tests; probe re-run |
| `5b120ec` | D6 + S7 local BLAST by database prefix, with provenance | external tool execution (bundled database) |
| `e1a48c6` | D7 stale BLAST artifact labeled as legacy (file kept) | docs |
| `7a4a718` | Handover notes (working practice) | docs |
| `3514617` | D9 CRISPR on-target site not an off-target | internal tests |
| `ca2ec39` | D8 CRISPR cut site per system, cited | internal tests; PubMed (PMIDs 22745249, 25830891, 30166441, 26422227) |
| `8145556` | D10 CRISPR config fields act; D11 exports in `exports/crispr` | internal tests (defaults unchanged) |
| `e47cb85` | D19 pipeline `run_orf` search settings; fixed `completed` shape | internal tests (equal to `core.api.find_orfs`) |
| `f3a09c1` | D15 planner entry labeled "not implemented" | end-to-end `main.py` |
| `c71963b` | S15 variant→mutation bridge works; S16 flag leak | internal tests; pyright −2 |
| `c743beb` | S17 cloning digest uses the top-strand cut | agreement with `restriction.perform_digest` for 231 pairs, linear and circular |
| `6c20e69` | S19 start-codon classification and scores use the genetic code | internal tests; golden unchanged |

## 4. Findings status

**Corrected:** D1–D15 (D14 with all sub-findings F1–F12), D19; S7, S15,
S16, S17, S19, S20.

**Open – defect work (to resume first, in this order):**

1. **D16** primer duplicates and the live `_v5`–`_v11` ranking chain.
   A safe first step is provable by bytecode: rename the first
   `calculate_tm_nearest_neighbor` definition to the private name the
   adapter already captures (`_calculate_tm_nearest_neighbor_v3_raw`)
   and drop the alias, since module-level lookups resolve at call time.
   Flattening the v5–v11 chain changes executed code and needs the
   independent review named in `docs/AUDIT_STATUS.md` (open validation
   item 7); do not do it without that.
2. **D17** (shared terminal visual language, `[0]` semantics) is the
   core of **Stage F** (shared renderer + `docs/TERMINAL_UI_STANDARD.md`).
   Decide with the owner whether D17 is done inside Stage F.
3. **D18** ("confidence" vocabulary in data and exports) changes the
   golden snapshots and is the core of **Stage G**; it needs a documented
   golden-master change. Same decision as for D17.
4. Remaining S-findings: review S1–S6, S8–S14 (audit section 5) and
   correct the confirmed ones; S18 (planner runtime, performance, P2-6);
   S21 (GenBank ORIGIN keeps only letters); S22 (variant reference
   alphabet unchecked). **S23, S3, S4 and S5 are done** – see section 7.
5. Not defects but recorded limitations: RNA module's three-frame codon
   view uses table 1 (P1-9); RID resume (P1-7); planner P3-1; deletion of
   `blast_results/TEST_8_REFERENCE_blast.json` is the owner's decision.

**Then:** Stages F–H and Stage I as specified by the Master Prompt
(`TERMINAL_UI_STANDARD.md`, `STYLE_GUIDE.md`, README, final report A–K).

## 5. Behavior changes a user will notice

- A DNA record containing IUPAC ambiguity codes other than N (e.g.
  `TEST_8_AMBIGUOUS` in `examples/dna.fasta`) asks once for its type
  (DNA/RNA/protein); `[moltype=genomic DNA]` in the header declares it.
- Protein FASTA/GenPept records open a protein menu (analysis, BLAST,
  report) instead of the DNA overview; the record list shows `aa`.
- RNA records are analyzed in the DNA alphabet with a printed note.
- `core.api` refuses wrong or undetermined sequence types
  (`SequenceTypeError`); every entry point accepts `sequence_type=`.
- CRISPR guides show the modeled cut; relative CRISPR exports go to
  `exports/crispr/`; the planner entry says "not implemented".

## 6. Continuation steps

1. Read `docs/PHASE6_PROGRESS.md` (handover notes: ReadOnly files, line
   endings, `orf.py` bare-LF lines, pyright comparison) and this file.
2. Verify: `git status` clean at `6f7987b`; run the suite (expect 228 OK);
   pyright 132 errors.
3. Resume with item 4.1 (D16, safe step only), then the remaining
   S-findings of item 4.4 (S1–S6, S8–S14, S18, S21, S22), then agree on
   4.2/4.3 (D17/D18) with the owner before Stages F–G.
4. One commit per logical package; restore ReadOnly after each commit.
   Attribute each commit to the model that actually performed the work.

Tooling hint for this environment: the Bash tool collapses `\\` in
heredocs; write edit scripts containing backslashes to a file (Write
tool) and run them, and write files with `newline="\n"` on Windows to
keep LF files LF.

---

## 7. Addendum: restricted CRISPR assignment (2026-10-31, Opus 5)

A model safeguard stopped the Opus 5.5 session at the start of the CRISPR
investigation covering S23, S3, S5 and S4. The owner assigned that single
blocked task to Opus 5, with an explicit instruction not to expand into
other Phase 6 work. This addendum records what that assignment produced.

**Scope actually permitted and completed.** All four findings proved
investigable and correctable without any restricted operation. No
safeguard was encountered during the work, nothing was bypassed, and no
task was reframed or split to work around one. The scope was limited by
the assignment, not by a capability limit.

| Finding | Outcome |
|---|---|
| S23 | Confirmed by execution and corrected. `seed_sequence` stored the guide's 3' end for every system; for Cas12a (5' PAM) that is the PAM-distal end. The new shared `pam_proximal_seed()` now serves both the stored seed and the local specificity screen, which had already been correct and therefore disagreed with the stored value. |
| S3 | Confirmed and corrected. `edit_in_rtt` only compared RTT lengths; it is now a coordinate check that the RTT window covers the edit site. The emission filter already checked the length separately, so the change tightens it without loosening anything. |
| S4 | Resolved by documentation. No threshold value was changed; each is now either cited or explicitly labeled a Synapse convention. |
| S5 | Partly confirmed. The "position-independent" concern is **not** a defect: Pol III terminates wherever it meets the T-tract. The genuine gap was the missing name and source, now supplied. The threshold was deliberately left unchanged. |

**Evidence.** Literature verified through PubMed rather than from memory:
Sternberg 2014 (PMID 24476820), Komor 2016 (PMID 27096365), Gaudelli 2017
(PMID 29160308), Anzalone 2019 (PMID 31634902), Verosloff 2021
(PMID 33811918). Nine regression tests were added. No external tool
execution was required for these findings, and none was claimed.

**Remaining limitations, deliberately not addressed here.**

1. The prime-editing RTT is still the **reference** template and does not
   encode the alternative allele, so it is not a ready-to-use pegRNA RTT.
   This is now stated in every candidate's notes and in METHODS. Emitting
   an edited RTT is a functional extension beyond S3 and needs an owner
   decision, because it changes what the module offers.
2. The exported field is still named `edit_in_rtt` although it now means
   "the RTT window spans the edit site". Renaming an export key belongs to
   the D18 terminology work package, which was out of scope here.
3. The seed length (8 nt) and the GC/homopolymer thresholds remain Synapse
   conventions. Deriving tool-specific values would require a benchmark
   that was not part of this assignment.
4. S5 was not turned into a position-dependent or promoter-aware check,
   because the evidence does not support a sharper rule: termination
   efficiency at U6 depends on tract length and on upstream RNA secondary
   structure, neither of which Synapse models.

**Incident during this assignment.** While updating this report, a
one-line Python command truncated it to an empty file: the destination was
opened for writing, which truncates, before the argument that read it was
evaluated. The file was unmodified in Git at that moment, so
`git checkout -- docs/PHASE6_HANDOVER.md` restored it completely and no
committed work was lost. Recorded here because the master prompt requires
that problems be reported rather than quietly repaired.

**State at the end of this assignment.** HEAD as recorded in section 1,
working tree clean, 228/228 tests pass, pyright 132 errors with no new
diagnostics, golden master unchanged, ReadOnly restored on
`analyses/crispr.py`.

**Next.** The main Phase 6 programme remains with Opus 5.5 and resumes at
section 6 above, with S23, S3, S4 and S5 removed from the open list.

---

## 8. Addendum: second Opus 5 assignment (2026-10-31)

The Opus 5.5 session was stopped by a safeguard again when it tried to
resume the programme, and the owner assigned the continuation to Opus 5.
This addendum records that block of work. It ends at a deliberate handover
point, not because the programme is finished.

**A note on how this work was split.** Two Opus 5.5 responses were stopped
by a classifier and the work was then routed to Opus 5. The tasks
themselves are ordinary scientific software maintenance, and no safeguard
was encountered while doing them. The stops look like false positives on
CRISPR vocabulary in a bioinformatics codebase. That is worth reporting to
Anthropic rather than worth designing around, and this session did not
sort or schedule work by expected classifier behavior.

### Completed in this assignment

| Commit | Work | Evidence |
|---|---|---|
| `c52f9bb` | **D16 part 1.** The Tm model is defined under its private name instead of carrying the public name and being shadowed by the compatibility adapter, so no module-level name is defined twice. | SHA-256 fingerprint of the compiled bytecode of all 129 module-level callables: no name added or removed, no entry changed. pyright's `reportRedeclaration` gone. |
| `5a3220e` | **S21, S10, S8, S22.** Four paths that discarded information or swallowed a failure without saying so now report it: symbols dropped from GenBank ORIGIN, feature locations the reader cannot parse, a failing `to_dict()` during report serialization, and a variant reference that is not DNA. | 7 regression tests; the bundled GenBank fixture produces no new notes. |
| `8c47b17` | **S1, S2, S9, S11, S12.** Two suspected findings assessed as not defects and locked by tests, one quantified by measurement, two resolved in documentation. | Execution on constructed minus-strand cases; `tracemalloc` at three sizes. |

### Findings resolved, with the honest outcome

- **S1, S2: not defects.** The minus-strand Kozak frame is the correct
  frame and reports no coordinate; nested minus-strand ORFs collapse
  correctly. Both are now covered by tests so the behavior cannot drift.
- **S9: confirmed, algorithm unchanged.** Measured at ~45 bytes per
  alignment matrix cell, so ~180 MB at 2,000 nt and ~4.5 GB at 10,000 nt.
  The note states these figures. A banded or linear-space alignment is
  still open as P2-6.
- **S11: still unverified, by choice.** Live literature-provider behavior
  needs external requests, which require the owner's approval under the
  resource policy. Offline fixtures continue to cover parsing, dedup and
  retries.
- **S12: resolved.** pyright runs against the baseline for every package.

### Still open, as independent work packages

1. **D16 part 2** - flatten the live `_v5`-`_v11` primer ranking chain.
   Blocked on purpose: it changes executed code and needs the independent
   review recorded as open validation item 7 in `AUDIT_STATUS.md`.
2. **S13** - `main.py` has module-level control flow and no
   `if __name__ == "__main__":` guard, so importing it starts the
   application. Fixing it means indenting roughly 700 lines into a
   function. Mechanical but large, with real regression risk for the smoke
   tests; it belongs with Stage E and deserves its own package.
3. **S14 / P2-2** - `core/api.py` still lacks cloning, variant filtering,
   mutation reporting, BLAST and literature entry points, and the analysis
   modules still own their menus.
4. **S6 / P1-7** - resuming an NCBI BLAST run from a stored RID. The
   primitives exist; the feature and its menu entry do not.
5. **S18** - cloning-planner runtime.
6. **P1-9** - the RNA module's three-frame codon view still uses table 1.
7. **Stages F, G, H, I** of the Master Prompt: terminal UI standardization
   (with D17), scientific output standardization (with D18, which changes
   the golden master and needs a documented justification), the style
   guide, editorial review and README, and the final verification with the
   A-K completion report.

### State at the end of this assignment

HEAD `8c47b17`, working tree clean, 240/240 tests pass, pyright 131 against
the 134 baseline with no new diagnostics, golden master unchanged, ReadOnly
restored on all 15 protected files.

Phase 6 is **not** complete. The defect and S-finding track is close to
done; stages F through I have not been started.

---

## 9. Addendum: third Opus 5 assignment (2026-10-31)

Six assigned work packages, each its own commit. Phase 6 is **not**
complete: stages F to I have not been started.

### Commits

| Commit | Package | Evidence |
|---|---|---|
| `dfe773d` | **D16 part 2.** The superseded primer generations are gone. | Investigation first, as the audit required: AST pass over every module-level binding in source order, project-wide reference scan, `sys.settrace` on a real design run. Bytecode fingerprint of all module-level callables: 129 -> 105, **nothing added, nothing changed**. 1,326 lines removed; pyright 134 -> 112. |
| `8d28586` | **P1-9.** The RNA codon view uses the selected ORF's genetic code. | For a table-2 ORF the +1 frame preview now equals `translate_orf_sequence(orf)` exactly. |
| `87c3fb0` | **S18.** Recognition-site scans are cached. | cProfile: 5,002 scans for 16 distinct combinations. Measured 11.8 s -> 0.24 s at 8 enzymes; 10 enzymes, timed at 100-110 s in the audit, 0.57 s. |
| `6a8c00c` | **S14.** The API covers cloning, variant filtering, mutation reports, BLAST and literature. | 16 entry points. The two that make external requests refuse to run without an explicit confirmation argument. |
| `ad843c5` | **S6 / P1-7.** An NCBI search can be resumed from its request ID. | Tested against a stubbed client that asserts `submit()` is never called, so no request left the machine. |
| `7008cd8` | **S13.** The session runs from `main()` behind a `__main__` guard. | `git diff -w` reduces the whole change to the 12 added lines. |

### Decisions worth knowing

- **D16's investigation corrected the audit.** The audit assumed
  `design_primers_professional` was the v7 alias. It is not: a plain
  definition near the end of the module overrides it. Only that pass over
  binding order in source order revealed it, and it changed which code was
  identified as dead.
- **The API must not quietly drop the privacy gate.** The BLAST menu asks
  before sending a query to NCBI. An API cannot ask, so `blast_search`
  refuses a remote engine without `send_to_ncbi=True`, `literature_search`
  requires `send_to_providers=True`, and `blast_resume` requires
  `contact_ncbi=True`. A script cannot transmit sequences by accident.
- **Nothing is invented for a resumed run.** BLAST output does not contain
  the query sequence, so a resumed run keeps an empty sequence and the
  reported length; when the XML states no length, the query says so.
- **S18 changed no science.** The cached scan is a pure function of
  sequence, recognition sequence and topology, and callers get fresh dicts
  so none can corrupt the cache.

### Still open

1. **Stage F** - terminal UI standardization, with **D17** and the
   remainder of **P2-2** (the analysis modules still own their menus). The
   BLAST menu's new entry is numbered `[14]` inside the SEARCH block, so a
   coherent renumbering belongs to this stage.
2. **Stage G** - scientific output standardization, with **D18**. This
   changes the golden master and needs a documented justification.
3. **Stage H** - `docs/STYLE_GUIDE.md`, editorial review, README.
4. **Stage I** - final verification and the A-K completion report.
5. **S11** - live literature-provider behavior, pending the owner's
   approval for external requests.
6. **S24 (new)** - BLAST XML is parsed with `xml.etree.ElementTree`.
   Modern CPython does not resolve external entities, so XXE is not
   reachable, but `defusedxml` would be a new runtime dependency needing
   its own justification. Recorded, not changed.
7. Deferred audit items unchanged by this assignment: P0-4 (systematic
   REBASE re-check), P2-3 (CFD scoring - still rejected without a
   published matrix), P2-4, P2-5, P2-6 (the mutation alignment is still
   quadratic; S18 addressed only the planner), P2-7 (Type IIS), and the
   P1-9 sub-items organism mode, 5'-partial ORFs and ORFs across the
   origin.

### Verification at the end of this assignment

Full suite 267/267, 0 skipped. Golden master 12/12 and the built-in module
self-tests pass separately. pyright 112 against the 134 baseline with no
new diagnostics. A cross-module run (DNA -> ORF with table 11 -> RNA ->
protein -> report) confirmed the ORF's table reaches the RNA view and the
+1 preview matches the protein translation. Working tree clean, ReadOnly
restored on all 15 protected files.

---

## 10. Addendum: Stage F complete (2026-11-01)

Stage F is done and verified. Stages G, H and I have not been started.

### Commits

| Commit | Work |
|---|---|
| `5db21ab` | One separator width through the new `ui/style.py`. The inventory counted **11 widths** across 389 literal sites; 382 literals in 16 files now call `major_rule()`/`minor_rule()` and the inventory reports **zero** remaining distinct widths. |
| `b1fa927` | `[0]` means one thing - leave the screen - as Back, Cancel or Exit (**7 meanings → 3**). Message prefixes unified to `ERROR:`/`WARNING:`/`NOTE:`/`OK:`. `docs/TERMINAL_UI_STANDARD.md` written; `tests/test_terminal_ui.py` enforces the standard project-wide. |
| `10431d6` | The BLAST menu renumbered: its SEARCH block read 1, 2, 3, 14 after the resume entry was added, and now reads 1-4, with the whole menu running 1-14. |
| `f803577` | S24 measured rather than assumed, and every outstanding audit item classified against the Master Prompt. |

### What Stage F deliberately did not do

- **No wording of a scientific statement changed.** Stage F moved frames,
  keys and prefixes only. Terminology belongs to Stages G and H.
- **Menus stay in the analysis modules.** `core/api.py` already covers
  every analysis, so nothing is forced through the menus; relocating them
  is the rest of P2-2 and a structural change of its own.
- **One deviation is recorded, not hidden.**
  `mutation._prompt_sequence` offers no `[0]`, because a real cancel needs
  a cancel contract and four changed call sites - control flow, not
  presentation.

### S24, measured

Not a clean pass and not a clean fail. **XXE is not reachable**:
`xml.etree.ElementTree` does not resolve external entities. **Entity
expansion is reachable but capped** by Expat's amplification guard (1 MB
expands in 0.02 s; 10 MB is blocked in 0.09 s). **The cap is
version-dependent** - an Expat 2.4.0 feature - which matters because
Synapse targets Python 3.10. The test skips rather than falsely passing on
an older build. `defusedxml` remains unjustified; the dependency-free
hardening, if ever needed, is to set the limits through
`xml.parsers.expat`.

### Remaining work

1. **Stage G** - scientific output standardization, with **D18**: the
   "confidence" vocabulary still in the data model and exports. This
   changes the golden master and needs a documented, scientifically
   justified snapshot update. It is a self-contained package and was not
   started here rather than left half-applied.
2. **Stage H** - `docs/STYLE_GUIDE.md`, editorial review, README.
3. **Stage I** - final verification and the A-K completion report.
4. The classified audit items in `docs/PHASE6_PROGRESS.md`: **P0-4** and
   **S11** are blocked on the owner's approval for external access;
   **P2-3**'s CFD scoring stays rejected without a published matrix;
   **P2-4**, **P2-6** (mutation alignment), **P2-7**, **P2-2** and the
   **P1-9** sub-items are open and unblocked.

### Verification

Suite 283/283, 0 skipped. Golden master and the built-in module
self-tests pass separately. The 12 interface-standard checks pass, and
were verified against deliberate violations before being relied on.
pyright 112 against the 134 baseline, no new diagnostics. Working tree
clean, ReadOnly restored on exactly the original 15 files.

## 11. Addendum: Stages G and H (2026-11-02, Opus 5.5)

Resumed from `aa5703b` on the owner's instruction to complete Stage G
(with D18), Stage H and Stage I. No external request was made; S11 stays
pending the owner's approval.

### Commits (oldest first, after `aa5703b`)

| Commit | Package |
|---|---|
| `317097d` | D18: heuristic scores renamed from "confidence" in data, exports and messages; ORF golden key renamed, values proven identical |
| `bc89394` | D18 part 2: the sequence-type classifier reports an `evidence` level |
| `ebf4d52` | One notation per unit (%, °C, cal/(mol·K), Da); primer positions shown and read 1-based |
| `cef464a` | A result's own warnings reach the exported summary (6.4) |
| `81ccc55` | Restriction parameters and the enzyme data version recorded (6.4) |
| `1a8f2bc` | Editorial review part 1: American spelling, one wording per message, tone |
| `4e5dbdc` | BLAST E-values with three significant figures everywhere |
| `97086ad` | `docs/STYLE_GUIDE.md` with a checked enforcement table |
| `ccf4ccb` | README rewritten to 15.1 |
| `1269335` | "Sequence complexity" defined with its limit; METHODS section for the sequence overview |

### Golden master

One intentional change (`317097d`): in `orf.json` and
`orf_legacy_baseline.json` the key `confidence` became `heuristic_score`
(6 and 23 occurrences). With the key mapped back the files equal the
previous ones, and the unmodified golden test confirmed the code
produces them. Hashes before are in `docs/CHANGES.md`. No other snapshot
changed.

### Recorded rather than changed

- `site_position` is 1-based in the restriction map and 0-based in digest
  records; it is stored in the golden master and the pUC19 comparison.
- Menu-path sequence-overview and cloning results record no input
  parameters.
- HGVS parentheses for predicted protein consequences are not used.
- The "complexity" field names are kept; the meaning is now documented.

### State at this checkpoint

Suite 298/298, 0 skipped; pyright 112 against the 134 baseline, no new
diagnostics; working tree clean; ReadOnly on exactly the original 15
files.

### Remaining work

Stage I: full verification, integration workflows (10.x), the local
reference comparisons that the environment allows (primer3-py 2.4.0 and
Biopython 1.88 are installed; `METHODS.md` still says Primer3 could not
be installed), dependency and compatibility review, documentation
consistency (`AUDIT_STATUS.md` is dated 2026-10-30 and its test count and
P2-6 status are stale), the feature matrix, git verification and the
A-K completion report.
