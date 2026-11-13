# Phase 6 completion run – continuation record

**Final state (2026-11-13).** The completion run ended with the final
verification, the independent audit and the new final report
(`docs/PHASE6_FINAL_REPORT.md`). Verdict: Phase 6 is PARTIALLY COMPLETE;
every safely executable requirement is done and verified, and the open
items are listed there in section J with what each needs (owner approval,
system-wide tools, large data, a licensed data file, or a specification).
The sections below are the history of the run.

Session of 2026-11-11, stopped by the owner after work packages A and B so
that the next session starts with full context. **Phase 6 is not
complete.** Start the next session from this file.

State at stop: `main` at the commit that adds this file, working tree
clean, ReadOnly restored on the 15 files, suite 421 tests OK, golden master
unchanged, pyright 111 errors with no new diagnostic against the baseline
of `ab4bd2a` (134).

## Done in this run

| Commit | Package | Item |
|---|---|---|
| `68e59d1` | A | Literature and BLAST clients close the response an HTTPError carries (the "ResourceWarning in a retry test" was a client leak) |
| `8379ac2` | A | Mutation sequence prompt gets `[0] Cancel` (`SequencePromptCancelled`); ten messages drawn as fields become ERROR/NOTE lines, with a rule; `MissingFileTests` (11.4) |
| `b4b16bb` | A/B | Cloning planner search cache (P2-6): 21.7 s -> 3.2 s, results byte-identical |
| `4939688` | B | Type IIS enzymes BsaI, BsmBI, BbsI, SapI on both strands in digests (P2-7); 140 cuts vs Biopython 1.88, 0 differences |
| `b7eb1de` | B | HGVS `p.(...)` predicted notation, used in the mutation report and variant `protein_change` |
| `1ae2833` | B | One coordinate convention for restriction site and cut records |
| `fedb750` | B | G·U wobble count in the RNA hairpin screen as a separate figure (P2-5) |
| `f020a1d` | B | ORFs across the origin of circular sequences and 5'-partial ORFs (P1-9) |
| `27255df` | B | Primer hairpin/dimer screens compared with primer3-py (rho -0.78 / -0.74 / -0.51) |
| `2c68e61` | B | Primers and guides checked against variants, `core.api.cross_check_variants` (P2-4) |
| `f9d4356` | B | Mutation alignment: one traceback byte per cell (46 -> 1.1 bytes) |
| `e5bb456` | B | Exact banded alignment: 10 kb similar sequences 23 s -> 0.18 s |
| `15e2ff1` | B | Banded Levenshtein similarity (10 kb 0.10 s); hairpin screen 135.6 s -> 25.6 s at 10 kb primer design |

Each commit message holds the evidence (before/after measurements,
byte-identical comparisons, tests that fail on the previous code).

## Done on 2026-11-12 (second run)

| Commit | Item |
|---|---|
| `bf4d096` | Approval policy of the Master Prompt adopted (owner decision); handover spelling fixed (the suite failed at `bd91a24`, 420/421) |
| `aba3fb0` | T1: primer design 10 kb 71 s -> 10 s CPU time (same machine and load), results identical (tier screen with an exact upper bound, faster stem/hairpin/repeat screens, caches) |
| `c6f0df9` | D20: the V5 structural penalty read missing fields and was always 0; D21: tm_method "nearest-neighbor" crashed candidate generation |
| `99150da` | T2: primer design around a chosen restriction site (10.2) in primer.py, `core.api`, the primer menu (entry 2; Single Primer Analysis is now 3) and the pipeline parameters |
| `50dea66` | T3/P2-2: the dialogs of eleven analysis modules moved to `ui/<module>_menu.py` (55 screens byte-identical); `core/textformat.py`; `SynapsePipeline.select_orf(chooser=...)`; no analysis module imports `ui` |
| `b1a012f` | T4: the sequence overview and the three cloning dialogs record their parameters (mode, enzymes, topologies, thresholds) in the pipeline and the report |
| `757958f` | T6: pyright 108 -> 58; D22 (orf_reference_cds AttributeError), D23 (pair risk assessment missing from exports) |
| `ae87adb` | T7: coverage 68 % -> 70.9 % (`docs/TEST_COVERAGE.md`); D24 (BLAST text export and query file written with CR CR LF) |
| `bf4b9ee` | T9: `docs/TOOL_AND_DATA_LEDGER.md`; suite on Python 3.10.11 (452 OK); pUC19 L09137.2 fixture with Biopython comparison in the suite |
| `dabc2f1` | T11: BLAST hit selector accepts only listed numbers; query-type question re-asks |
| (this commit) | README brought up to date; live literature check recorded in METHODS; `max_results` per provider documented in `core.api` |

## Resumed on 2026-11-13 and finished

| Commit | Item |
|---|---|
| `a5e5a94` | D25: BLAST polling survives brief network failures and names the RID; primer3-py comparisons skip when Windows application control blocks the package |
| `cdbf69c` | T12: variant check in the primer and CRISPR menus |
| `f58d866` | T8: CHANGES, FEATURE_STATUS, AUDIT_STATUS, progress record |
| final commits | Last padded labels; Ruff pyflakes pass (40 unused imports removed); REBASE 611 comparison 26/26; tool evaluation (7.4); final report A-K with requirement matrix and independent audit |

## Paused here on 2026-11-12 (owner asked to stop, historical)

- **D25, found, not yet fixed:** a live NCBI BLAST search (public pUC19
  fragment 401-700, database `nt`) was submitted and got an RID, then one
  status poll failed with a transient DNS error (`getaddrinfo failed`) and
  the whole search aborted. `NCBIBlastClient.await_ready` does not retry
  connection errors while polling, and the error message does not contain
  the RID, so the user cannot resume a search that is still running at
  NCBI. Fix: retry transient connection errors in `await_ready` (bounded,
  with back-off) and name the RID and the resume entry in the final error;
  test with a client whose `_request` fails once. Then repeat the live run
  (script idea: read `tests/fixtures/pUC19_L09137.2.fasta` directly, do not
  import `tests._support`, which blocks the network on purpose).
- Remote BLAST does not offer NCBI's `core_nt` (NCBI's default since 2024);
  decide after the live run whether to add it.
- P2-1 decision: keep the v5-v11 names (one live implementation per role;
  renaming changes no behavior). Record as PARTIALLY COMPLETE with reason.
- Still open: T12 (P2-4 in the menus), T13 (document), T8 documentation
  reconciliation for all commits since 2026-11-03 (CHANGES, FEATURE_STATUS
  incl. new file locations ui/*_menu.py, AUDIT_STATUS, PHASE6_PROGRESS,
  FINAL_REPORT J), T10 final verification, independent audit and final
  report A-K.

Measurement note: the owner's machine ran other programs (CPU load 74-87 %)
during this run; compare timings only within one run.

## Open, in the planned order

### B leftovers (scientific / performance)

1. **Primer design at 10 kb: 25.6 s**, target under 10 s (P2-6 "primer
   pre-filter"). Largest cost now `_alignment_stem` (self- and heterodimer,
   `analyses/primer.py` ~line 380), then `normalize_dna_sequence` called
   ~2 million times. Approach as before: snapshot `design_primers_professional`
   on 600/900/1,200 bp (serialized dataclasses), optimize, prove identical.
2. **Restriction-site-driven primer design (10.2 workflow)**: no wiring from
   a chosen restriction site to the primer target region.
3. **Merging the surviving primer generations (P2-1)**: deliberately not done
   so far; decide whether it is required.

### C – architecture

4. **P2-2: menus out of the analysis modules.** All twelve analysis modules
   import `ui.style` because their menus live in them (133 top-level
   definitions, ~8,900 lines). Plan: move the style-using functions to
   `ui/` modules (AST tool that collects free names and imports them from
   the analysis module), update `main.py`, tests (`render_first_screen`,
   `test_rendered_menus`, `INTERFACE_SOURCES`), then a test that imports every
   analysis module with `ui` refused. Exceptions to settle:
   `mutation.build_mutation_report` (uses `field_line` for report text),
   `orf.select_orf`, the self-test printers.
5. **Pipeline menu-path parameters** for sequence overview and cloning.

### D – quality assurance

6. pyright triage (111 errors, reachable `None` paths first).
7. Line coverage with the standard-library `trace` module (`coverage` is not
   installed; installing needs approval).
8. Tool and data ledger (Master Prompt 1.3.1, 11.8): Biopython 1.88,
   primer3-py 2.4.0, BLAST+ 2.17.0, pyright 1.1.414, the pUC19 and other
   fixtures with source, date, licence, checksum.

### Documentation not yet updated for this run

`CHANGES.md` (one entry per commit above),
`FEATURE_STATUS.md` (cloning planner, Single-Enzyme, Type IIS, ORF partial
and circular rows 6.5 and "Circular / origin-spanning ORFs", RNA hairpins,
mutation report notation and alignment, variant bridge, new API function,
test count), `AUDIT_STATUS.md` (P1-9, P2-4, P2-5, P2-6, P2-7),
`PHASE6_FINAL_REPORT.md` section J (remove the rows done above) and the
README capability list (Type IIS, circular ORFs). `METHODS.md` and
`TERMINAL_UI_STANDARD.md` were updated in the commits.

### Final steps

9. Full verification (suite, golden master, self-tests, rendered menus,
   architecture checks, pyright, diff/whitespace).
10. Independent final audit against the whole Master Prompt and the final
    report A-K with the COMPLETE / PARTIALLY COMPLETE / UNVERIFIED / BLOCKED
    matrix.

## Approval policy since 2026-11-12

The owner decided that the Master Prompt's approval rules apply (§1.3,
§14.5); see `docs/PROJECT_SYNAPSE_COMPLETE_HANDOVER.md`, Part X. Project-local
installs, small downloads and consulting public sources are allowed without
asking. Approval is needed for very large downloads, system-wide changes,
costs, uploading private data or sequences, and irreversible operations.
Several rows below were therefore re-assessed in the handover, Part XIII.

## Blocked (needs the owner, as recorded on 2026-11-11)

| Item | Needs |
|---|---|
| Python 3.10 runtime test | installing a 3.10 interpreter (only 3.14.6 present) |
| `bcftools norm`, ORFfinder comparisons | WSL or a Linux host (WSL not installed) |
| Live NCBI BLAST run, RID resume, live literature providers (S11) | approval for external requests |
| REBASE direct re-check, NEBcutter, PubMed check of four METHODS references, Guruprasad 1990 table | approval for external requests |
| CFD off-target scoring (P2-3) | the published matrix as a data file with citation |
| `coverage` line measurement | approval to install (standard-library `trace` is the fallback) |
| Licence, `blast_results/TEST_8_REFERENCE_blast.json` | owner decisions |

## Working notes

- Helper scripts lived in the session scratchpad and are not in the
  repository; `patchlib.apply()` (exact, CRLF-preserving replacements),
  `eol.py`, `pyright_diff.py`. The pyright baseline is reproducible:
  check out `ab4bd2a` in a worktree, point `venvPath` at this project and
  write `file|rule|message` lines of `python -m pyright --pythonversion 3.10
  --outputjson` (134 errors).
- Run pyright as `python -m pyright` (the `.exe` is blocked by application
  control). Clear and restore ReadOnly with `Set-ItemProperty ... IsReadOnly`.
