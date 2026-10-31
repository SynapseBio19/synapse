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
| HEAD | `6c20e69` (plus the commit that adds this file) |
| Working tree | clean |
| Tests | 219/219 pass, 0 skipped (`python -m unittest discover -s tests -t .`, ~25 s) |
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
   alphabet unchecked); **S23** (Cas12a `seed_sequence` is the PAM-distal
   end; the off-target screen already uses the right end – a one-line
   fix in `_candidate_from_oriented_window` plus a test).
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
2. Verify: `git status` clean at the commit that adds this file; run the
   suite (expect 219 OK); pyright 132 errors.
3. Resume with item 4.1 (D16, safe step only), then 4.4 (S23 and the
   remaining S-findings), then agree on 4.2/4.3 (D17/D18) with the owner
   before Stages F–G.
4. One commit per logical package; `Co-Authored-By: Claude Opus 5.5
   <noreply@anthropic.com>`; restore ReadOnly after each commit.

Tooling hint for this environment: the Bash tool collapses `\\` in
heredocs; write edit scripts containing backslashes to a file (Write
tool) and run them, and write files with `newline="\n"` on Windows to
keep LF files LF.
