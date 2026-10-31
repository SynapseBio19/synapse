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
| HEAD | `6f7987b` (see section 7; `642f31b` was the state at the original pause) |
| Working tree | clean |
| Tests | 228/228 pass, 0 skipped (`python -m unittest discover -s tests -t .`, ~28 s) |
| Golden master | unchanged since the audit; no snapshot was regenerated |
| pyright 1.1.414 (`--pythonversion 3.10`) | 132 errors; baseline at `ab4bd2a` was 134; no new diagnostic in any package; the two removed are the S15 call |
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
| `e1a48c6` | D7 stale BLAST artifact labelled as legacy (file kept) | docs |
| `7a4a718` | Handover notes (working practice) | docs |
| `3514617` | D9 CRISPR on-target site not an off-target | internal tests |
| `ca2ec39` | D8 CRISPR cut site per system, cited | internal tests; PubMed (PMIDs 22745249, 25830891, 30166441, 26422227) |
| `8145556` | D10 CRISPR config fields act; D11 exports in `exports/crispr` | internal tests (defaults unchanged) |
| `e47cb85` | D19 pipeline `run_orf` search settings; fixed `completed` shape | internal tests (equal to `core.api.find_orfs`) |
| `f3a09c1` | D15 planner entry labelled "not implemented" | end-to-end `main.py` |
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

## 5. Behaviour changes a user will notice

- A DNA record containing IUPAC ambiguity codes other than N (e.g.
  `TEST_8_AMBIGUOUS` in `examples/dna.fasta`) asks once for its type
  (DNA/RNA/protein); `[moltype=genomic DNA]` in the header declares it.
- Protein FASTA/GenPept records open a protein menu (analysis, BLAST,
  report) instead of the DNA overview; the record list shows `aa`.
- RNA records are analyzed in the DNA alphabet with a printed note.
- `core.api` refuses wrong or undetermined sequence types
  (`SequenceTypeError`); every entry point accepts `sequence_type=`.
- CRISPR guides show the modelled cut; relative CRISPR exports go to
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
| S4 | Resolved by documentation. No threshold value was changed; each is now either cited or explicitly labelled a Synapse convention. |
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
