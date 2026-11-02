# Phase 6 final report

Project Synapse, Phase 6 (2026-10-30 to 2026-11-02). Working copy
`Project Snyapse - Kopie`; the original project was not modified.
Governing specification: the integrated Phase 6 Master Prompt. Work by
Claude Opus 5 and Claude Opus 5.5. 58 of the 63 Phase 6 commits carry a
model attribution; the first five (`824ca5f`, `d4544ad`, `5c098e3`,
`47146b0`, `7334296`) carry none.

Corrected on 2026-11-02 after the closing review
(`docs/PHASE6_INDEPENDENT_REVIEW.md`): figures in sections A and G, the
evidence wording for the IG weight and the Python 3.10 finding, and
the open items in section J.

The detailed records behind this report are `docs/PHASE6_INITIAL_AUDIT.md`
(the findings), `docs/PHASE6_PROGRESS.md` (one row per work package, with
commit and evidence), `docs/CHANGES.md` (every change with before, after
and reason), `docs/FEATURE_STATUS.md` (every feature), `docs/METHODS.md`
(methods, defaults, validation) and `docs/PHASE6_HANDOVER.md` (the
execution history, including pauses and a restricted assignment).

---

## A. Executive summary

**Outcome.** The Phase 6 audit recorded 19 confirmed defects (D1-D19)
and, as the work went on, 24 suspected or unverified concerns (S1-S24).
All 19 defects are corrected, and every S item is assessed. Most
confirmed ones are corrected; S9 (alignment memory) is quantified and
documented but not corrected, S14 is resolved for the API while the
menus remain in the modules (P2-2), and S5 is addressed by naming and
sourcing the rule. Two were shown not to be defects, one (S24) is
measured with its residual risk documented, and one (S11) is blocked on
the owner's approval for external requests. Stages F (terminal
interface), G (scientific output), H (language and documentation) and
I (final verification) are complete. Three further defects that no
earlier test had caught were found and corrected (section B, items
25-27): the primer coordinates in Stage G, the IG weight and the
Python 3.10 blocker in Stage I.

**Main improvements.**

- Scientific: initiator codons, ORF translation and start-codon scoring
  follow the genetic code; restriction, cloning and the circular map
  honor topology and ignore ambiguous sites; CRISPR guides report a
  modeled cut site with its source and no longer count their own site as
  an off-target; the instability index uses the reference dipeptide
  weights; heuristic scores are no longer called "confidence".
- Input: a central sequence-type layer (`core/alphabets.py`) classifies
  every record as DNA, RNA or protein with an evidence level and asks the
  user when the letters cannot decide, instead of loading a peptide as
  DNA.
- Architecture: one input contract for menu and API; `core/api.py`
  covers every analysis without keyboard input; one shared presentation
  module (`ui/style.py`); the superseded primer generations are gone.
- Interface and language: one rule width and one meaning for `[0]`
  across all screens; since the Stage F completion the ORF, CRISPR, BLAST,
  literature, primer and variant screens are drawn by the shared
  renderers (the other modules still print their own titles and fields);
  one message prefix per category, one
  notation per unit, 1-based positions everywhere a user sees them, and a
  style guide whose mechanical rules the test suite enforces.
- Evidence: 314 tests (109 at the end of the previous cycle), including
  external reference comparisons with Biopython 1.88 and primer3-py
  2.4.0, BLAST+ 2.17.0 executions, and cross-module workflow tests.

**Actual readiness.** Synapse is a coherent, transparent research and
teaching workbench whose calculations agree with established reference
implementations wherever a comparison was possible. It is **not**
experimentally validated, has not been run on a Python 3.10
interpreter, and several modules remain heuristic by design (guide
ranking, primer structure screens, signal-peptide and transmembrane
screens). It should be used where its stated limitations are acceptable
and its results checked with established tools before experimental
decisions.

**Important remaining risks.** No run on the minimum supported Python
(only a static check); no live test of the four literature services or
of a remote NCBI BLAST run; enzyme definitions checked against
Biopython's REBASE data but not systematically against REBASE itself;
the mutation alignment needs ~4.5 GB for two 10-kb sequences; the
analysis menus still live inside the analysis modules.

---

## B. Scientific improvements

Every row: previous behavior -> corrected behavior. The reasoning,
references and test names are in `docs/CHANGES.md` under the ID. Levels:
**R** external reference comparison, **X** external tool execution,
**K** known-answer test, **I** internal test only.

| # | ID | Module | Previous | Corrected | Evidence |
|---|---|---|---|---|---|
| 1 | D1 | ORF | Initiator Met for any GTG/TTG start, whatever the code | Met only for initiation codons of the ORF's genetic code | R (Biopython, tables 1/2/4/11) |
| 2 | D2 | ORF, protein, API | Translation fell back to table 1 | Each ORF is translated with its own table; a conflicting table is an error | R, K |
| 3 | D3 | Restriction | Map and sites ignored circular topology | Origin-spanning sites found in every site-based result | R (pUC19, 132/132) |
| 4 | D4 | Cloning | Ambiguous (POSSIBLE) sites could define cloning ends | One shared matcher; only CONFIRMED sites define ends | K, R |
| 5 | D5 | Cloning | No circular support for vector or insert | Topology of both molecules honored | K |
| 6 | D6, S7 | BLAST | Local database by path always failed; provenance empty | Path prefix scanned; title/date/count via `blastdbcmd -info` | X (BLAST+ 2.17.0) |
| 7 | D7 | BLAST | A scientifically wrong result file shipped unlabeled | Labeled as a legacy record (deletion is the owner's decision) | I |
| 8 | D8 | CRISPR | No cut site anywhere | Modeled cut per system, with citation, in model, display and exports | K (cited conventions) |
| 9 | D9 | CRISPR | The guide's own site counted as an off-target | On-target excluded from the local specificity screen | K |
| 10 | D10 | CRISPR | Five configuration fields were never read | Each field acts; defaults keep previous results | I |
| 11 | D12 | Variant | Contig-name mismatch skipped variants silently | Flagged; opt-in `chr` prefix matching | K |
| 12 | D13, D14 | Input, all modules | A peptide whose letters are all IUPAC nucleotide codes loaded as DNA; two `normalize_dna` contracts | Central classifier with evidence levels; declared type from GenBank/FASTA; the user declares undecidable records; one input contract | R (`Bio.Data.IUPACData`), X (blastp) |
| 13 | D14-F6 | Protein | B/Z/J/X/U/O lost or shifted | Every position kept; masses and ranges defined | R |
| 14 | D18 | ORF, protein, cloning, primer, classifier | Heuristic scores and type outcomes called "confidence" | `heuristic_score`, `evidence`; golden key renamed with values proven identical | I (golden master) |
| 15 | S3 | CRISPR prime editing | `edit_in_rtt` checked only the RTT length | Checks that the RTT contains the edit | K |
| 16 | S4, S5 | CRISPR | Thresholds uncited; poly-T flag unexplained | Sourced or labeled Synapse conventions; Pol III/U6 consideration named and cited | literature |
| 17 | S15, S16 | Variant | Variant-to-mutation bridge never worked; REF_MISMATCH leaked into input | Bridge works; flag kept in the result | K |
| 18 | S17 | Cloning | Digest cut at the site start | Cut at the enzyme's top-strand cut | K, R |
| 19 | S19 | ORF | Start classification and scores used a fixed ATG/GTG/TTG set | Initiation codons of the genetic code | K |
| 20 | S20, S22 | BLAST, variant | J/U/O refused; any reference alphabet accepted | J/U/O accepted for protein; the variant reference must be DNA | X, K |
| 21 | S21, S10, S8 | GenBank, report | ORIGIN, unsupported locations and failing serialization dropped data silently | Each loss is reported in a note or the log | K |
| 22 | S23 | CRISPR | Seed always taken at the guide's 3' end, which is PAM-distal for Cas12a (5' PAM) | PAM-proximal seed for 3'- and 5'-PAM systems | K |
| 23 | S9 | Mutation | Memory cost of the alignment unstated | Measured (~45 bytes per cell) and stated in the note | measurement |
| 24 | P2-5 | Sequence overview | "Complexity" not defined | Defined as normalized composition entropy with its limit (repeats score 1.0) | K |
| 25 | IG weight | Protein | Instability weight for IG was -7.49 | 1.0, the value in Biopython's table (`ProtParamData.DIWV`, which cites Guruprasad et al. 1990); the paper's table was not consulted | R (Biopython ProtParam, 31 proteins; one secondary source) |
| 26 | Py 3.10 | CRISPR | `f"{'#':>3}"` made the module unimportable on 3.10/3.11 | Literal header, same output | static 3.10 check; the import failure follows from the pre-3.12 grammar and was not observed on a 3.10 interpreter |
| 27 | Primer coordinates | Primer menu | Positions printed 0-based half-open unlabeled; target read half-open | 1-based inclusive in display and input; API unchanged | K |

Independent scientific validation (experiments, curated benchmarks)
exists for none of these; see section K.

---

## C. Functional completeness

`docs/FEATURE_STATUS.md` is the complete inventory: every module, menu,
submenu, function and capability with an implementation status, a test
status and a validation status, updated to the end of Phase 6.

- **Implemented and accounted for**: all eleven main-menu modules and the
  session functions (`[P]`, `[R]`, `[S]`). No original feature was lost.
- **Corrected**: the rows in section B.
- **Implemented in Phase 6**: CRISPR cut sites; CRISPR configuration
  fields; on-target exclusion; RID resume (`NCBIBlastClient.resume`,
  BLAST menu `[2]`, `core.api.blast_resume`); API entry points for
  cloning, variant filtering, mutation reports, BLAST and literature;
  protein-record route; sequence-type declaration.
- **Tested**: every module is execution-tested; the DNA-to-protein,
  restriction-and-primer, CRISPR, variant and reporting workflows are
  integration-tested (section H).
- **Independently compared**: genetic codes, translation, protein
  descriptors, enzyme definitions, digests, the circular map, end
  compatibility, Tm, local BLAST.
- **Not implemented, stated as such**: the experiment planner (the menu
  says so), 5'-partial and origin-spanning ORFs, Type IIS enzymes, CFD
  scoring, genome-wide off-target search, restriction-site-driven primer
  design, an organism mode for ORF search.

---

## D. Architecture

**Before Phase 6.** Analysis modules combined science and menus; each
module checked its own input with its own rules (two `normalize_dna`
functions with different contracts); no layer decided what kind of
sequence a record was; the primer module carried generations of
superseded definitions behind shadowing aliases; the API lacked several
analyses; each module drew its own screens.

**After Phase 6.**

| Component | Role |
|---|---|
| `core/alphabets.py` (new) | Alphabets, normalization and the DNA/RNA/protein classifier with evidence levels. One contract for every entry point. |
| `core/models.py` | `SequenceRecord` carries declared type, molecule type and the classification; `AnalysisResult` lifts a payload's warnings into its summary. |
| `core/api.py` | Every analysis without keyboard input, with the same type checks as the menus. |
| `core/pipeline.py` | Session state; `run_orf` takes the menu's search settings; ORF-dependent results are invalidated on reselection; restriction and cloning results record parameters and the enzyme data version. |
| `ui/style.py` (new) | All presentation: rules, headers, fields, prompts, message prefixes. Imports nothing scientific. |
| `analyses/restriction.py` | One recognition-site matcher shared with cloning; cached scans. |

**Removed duplication.** D16: the primer module's shadowed definitions
and superseded generations (1,326 lines in part 2, after a reachability
review with `sys.settrace`; every survivor bytecode-identical). D13: one
normalization contract. D17: 389 separator literals replaced by two
helpers.

**API changes.** New: `plan_cloning`, `filter_variants`,
`mutation_report`, `blast_search`, `blast_resume`, `literature_search`.
Renamed data keys (D18): `heuristic_score`, `heuristic_classification`,
`*_heuristic_score`, `evidence`, `type_evidence`,
`sequence_type_evidence`; `calculate_orf_confidence` remains as an alias.
Restriction and cloning results gain `enzyme_database_version`.

**Pipeline changes.** D19 parity with the menu; fixed status shape;
stale-result removal verified end to end.

**Modified files.** 50 files over Phase 6 (11,892 insertions, 3,516
deletions). Added: `core/alphabets.py`, `ui/style.py`, six test modules
(`test_phase6_tier1`, `test_phase6_tier2`, `test_phase6_sequence_types`,
`test_terminal_ui`, `test_external_reference`,
`test_integration_workflows`), `docs/FEATURE_STATUS.md`,
`docs/PHASE6_*.md`, `docs/TERMINAL_UI_STANDARD.md`, `docs/STYLE_GUIDE.md`
and `blast_results/README.md`.

---

## E. Terminal UI standardization

Measured before Stage F (`docs/TERMINAL_UI_STANDARD.md`, section 1): 11
separator widths across 389 sites; 7 meanings of `[0]`; three error
prefixes. After: one width (70) from `ui.style`; `[0]` = Back, Cancel or
Exit (main menu only); `ERROR:`/`WARNING:`/`NOTE:`/`OK:`; `Invalid
selection.` for every unknown key.

| Module | Before | After |
|---|---|---|
| Main menu, session | 60/68-column rules; `[0]` "Choose another sequence" in the type dialog | 70 columns, groups in work order, `[0] Exit`; the type dialog's `[0]` is `Cancel` with its effect named |
| Sequence overview | 70/45 rules | shared rules and fields |
| ORF | own frames; `[0]` "Apply and return" in the settings | shared frames; `[0] Back` (settings apply as entered); the comparison table rule grows with its columns (documented exception) |
| RNA, protein | own frames | shared frames, prompts and messages |
| Restriction, cloning | own frames | shared frames; `[0]` Back or Cancel |
| Primer | 0-based unlabeled positions | 1-based positions, labeled |
| Mutation, variant | `[0]` "Enter sequence manually" | `[M]` for manual entry |
| CRISPR | 78-column frames | 70 columns |
| BLAST | 90-column frames; SEARCH block 1, 2, 3, 14 | 70 columns; menu 1-14 in work order |
| Literature | 88-column frames | 70 columns |

Enforced by `tests/test_terminal_ui.py` (12 checks, verified against
deliberate violations). The Stage F re-examination
(`docs/STAGE_F_UI_REVIEW.md`) found that these checks covered rules,
`[0]` and prefixes only and that no module used the shared renderers;
the six modules named there were then converted and are checked as
rendered screens (`RenderedScreenTests`). **Remaining inconsistencies**:
ten interface files not yet converted (`TERMINAL_UI_STANDARD.md`
section 9); the mutation
module's sequence prompt offers no `[0]` (a cancel needs a control-flow
change in four callers); `core/pipeline.py` prints status from inside
`core/` through two function-local imports; the menus still live inside
the analysis modules (P2-2).

---

## F. Scientific output standardization

- **Result structures**: every stored result is an `AnalysisResult` with
  parameters, input SHA-256, warnings and limitations; warnings a payload
  carries now reach the summary on every path.
- **Units** (one notation each, enforced): `52.38%`, `60.12 °C`,
  `kcal/mol`, `cal/(mol·K)`, `Da`, `nt`/`bp`/`aa`.
- **Precision**: 2 decimals for percentages, Tm and thermodynamics; E-values
  3 significant figures everywhere (`n/a` when absent).
- **Coordinates**: users see 1-based inclusive positions; internal
  conventions are documented per structure. One naming inconsistency is
  recorded, not changed: `site_position` is 1-based in the restriction map
  and 0-based in digest records (both stored in the golden master).
- **Interpretation**: the measured / calculated / predicted / inferred
  distinction and the interpretation structure are in the style guide;
  heuristic values are named heuristic.
- **Exports**: JSON (complete) and Markdown (summary) with software
  version, input provenance, per-result parameters, warnings,
  limitations and a payload hash. Restriction and cloning record the
  enzyme data version. Menu-path sequence-overview and cloning results
  still record no input parameters (recorded).

---

## G. Language and editorial improvements

All meaningful changes are tabulated with original and revised wording,
reason and category in `docs/CHANGES.md` (Stage H sections, Stage G
D18). In summary:

- **Terminology**: "confidence" removed from data keys and messages;
  "evidence level" for the type classifier; "heuristic evidence" for the
  signal-peptide screen; the primer warning names what was observed
  ("Every reported primer pair is rated MARGINAL or POOR ...").
- **Spelling**: American English throughout (86 words in 8 documents,
  9 text sites in 4 source files and 6 in tests), matching the
  identifiers; enforced by a test.
- **Recurring messages**: one wording per situation (`Invalid
  selection.`, `No DNA sequence selected.`).
- **Precision**: the mutation frame warning names the reading frame and
  states the inference as a likelihood; the complexity interpretation
  states what the value measures and where it fails.
- **Tone**: no exclamations in messages.
- **Documents**: README rewritten to section 15.1; METHODS gained a
  sequence-overview section and four references and lost a
  self-contradiction; AUDIT_STATUS and FEATURE_STATUS brought to the end
  state; broken Markdown tables repaired.
- **Shared conventions**: message prefixes and helpers in `ui/style.py`;
  terminology, units, writing, messages and code rules in
  `docs/STYLE_GUIDE.md`, whose enforcement table is itself tested.

An inventory of every interface string found no German text and no claim
stronger than the evidence: each "prove" and "guarantee" is a negation.

---

## H. Testing

Final run at `db74504` plus the documentation commit that adds this
report, on the development machine:

| Measure | Value |
|---|---|
| Tests run | 314 |
| Passed | 314 |
| Failed | 0 |
| Skipped | 0 (BLAST+, Biopython and primer3-py all present) |
| Golden-master snapshots | 13, all matching; one intentional change in Phase 6 (D18, key rename, values proven identical) |
| Built-in module self-tests | all passing (run by the suite) |
| Integration workflows | 7 tests (10.1-10.4, 10.7), all passing |
| pyright 1.1.414, `--pythonversion 3.10` | 112 diagnostics against the Phase 6 baseline of 134; no new diagnostic in any package |

**Scientific reference comparisons** (external, not experimental):

| Comparison | Reference | Result |
|---|---|---|
| Melting temperature | primer3-py 2.4.0 | 2,916 comparisons over Mg²⁺ 0-3 mM, max 0.0032 °C |
| Enzyme definitions | Biopython 1.88 `Bio.Restriction` | 22/22 |
| Digests, two topologies, three rotations | Biopython 1.88 | 132/132 (6-kb test sequence; pUC19 in the audit, 132/132) |
| End compatibility | Biopython 1.88 | 484/484 pairs |
| Protein descriptors | Biopython 1.88 `ProtParam` | pI, GRAVY, aromaticity, extinction identical or within 1e-4; instability within 2e-4 after the IG fix; mass within 2e-5 relative |
| Genetic codes and translation | `Bio.Data.CodonTable` | tables 1, 2, 4, 11 identical |
| Alphabets | `Bio.Data.IUPACData` | identical letter sets |
| Local BLAST | BLAST+ 2.17.0 | real runs, provenance complete |

**Environment**: Windows 11; Python 3.14.6 (`.venv`); Biopython 1.88,
primer3-py 2.4.0 (test-only); NCBI BLAST+ 2.17.0; Expat 2.8.1; pyright
1.1.414. The suite blocks all network access. Not run: a Python 3.10
interpreter (not installed), live remote services.

One minor observation: a literature retry test leaves a simulated
`HTTPError` for the garbage collector, which prints a `ResourceWarning`;
it does not affect results.

---

## I. Git history

62 commits in Phase 6 before this report, from `824ca5f` (initial audit)
to `db74504`, each one logical package with tests and documentation; the
report and the final handover are committed after them.
No history was rewritten.

| Stage | Commits | Content |
|---|---|---|
| Audit | `824ca5f` | Initial audit and feature matrix |
| Tier 1 | `d4544ad` ... `ab4bd2a` | D12, D4, D3, D5, D1, D2 |
| Sequence types | `b8f5ea8` ... `d40d328` | D14, D13 |
| Tier 2 | `5b120ec` ... `6f7987b` | D6/S7, D7, D9, D8, D10/D11, D19, D15, S15-S17, S19, S23/S3-S5 |
| Architecture | `c52f9bb` ... `7008cd8` | D16, S8/S10/S21/S22, S1/S2/S9/S11/S12, P1-9, S18, S14, S6, S13 |
| Stage F | `5db21ab` ... `f803577` | D17, BLAST menu, S24 |
| Stage G | `317097d` ... `81ccc55`, `4e5dbdc` | D18, units, coordinates, exports |
| Stage H | `1a8f2bc` ... `1269335` | editorial review, style guide, README, complexity |
| Stage I | `2851e11` ... `db74504` | reference comparisons, IG fix, workflows, status documents, Python 3.10 |
| Records | `02d2899`, `7a4a718`, `642f31b`, `bff3c7d`, `edbc28f`, `692cb20`, `aa5703b`, `4c31602` | progress and handover |

Working tree: clean after the final commit; no uncommitted changes. The
Windows ReadOnly attribute is set on exactly the 15 files that carried it
at the start. Not committed by design: `exports/`, `logs/`,
`data/blast_db/`, caches and local skill files (`.gitignore`).

---

## J. Remaining work

| Issue | Module | Why unresolved | Needs | Next action | Kind |
|---|---|---|---|---|---|
| Live literature services (S11) | literature | external requests need approval | owner's approval | one controlled run per provider | environmental |
| Remote NCBI BLAST run and live RID resume | BLAST | external request with a query sequence | owner's approval | one confirmed run recording RID and provenance | environmental |
| Python 3.10 suite run | all | no 3.10 interpreter | installation approval | install 3.10, run the suite | environmental |
| REBASE re-check of all 22 enzymes (P0-4) | restriction | external data source | approval for REBASE access | compare every entry with REBASE | environmental |
| Instability table against Guruprasad et al. 1990 | protein | primary table not consulted | the paper | confirm IG and FY | scientific |
| `bcftools norm`, ORFfinder comparisons | variant, ORF | no native Windows build / no local binary | WSL or Linux host | run both on reference sets | environmental |
| Primer hairpin/dimer vs Primer3 | primer | heuristic by design, not compared | none (primer3-py is installed) | compare with `calc_hairpin`/`calc_heterodimer` | scientific |
| CFD off-target scoring (P2-3) | CRISPR | needs the published matrix as cited data | the matrix | add as a data file with citation | scientific |
| Primer/guide vs variant cross-checks (P2-4) | primer, CRISPR, variant | not implemented | none | implement | technical |
| Mutation alignment memory (P2-6) | mutation | algorithmic change | none | banded or linear-space alignment | technical |
| Type IIS enzymes (P2-7) | restriction | feature absent | none | add with both-strand search | technical |
| Menus inside analysis modules (P2-2) | all | structural change of its own | none | move menus to `ui/` | technical |
| Organism mode, 5'-partial, origin-spanning ORFs (P1-9) | ORF | not implemented | none | implement | technical |
| Experiment planner (P3-1), genome specificity (P3-2) | planner, primer, CRISPR | deferred | genome data for P3-2 | design first | technical |
| HGVS parentheses for predicted consequences | mutation | convention not applied | none | write `p.(R12C)` and test | scientific |
| `site_position` 1-based in map, 0-based in digests | restriction | changes stored values | none | rename with a documented snapshot change | technical |
| Menu-path parameters for sequence overview and cloning | pipeline | which inputs to record is a per-module choice | none | name them explicitly | technical |
| `[0]` in the mutation sequence prompt | mutation | needs a cancel contract in four callers | none | add and test | technical |
| Restriction-site-driven primer design | primer, restriction | no wiring | none | pass a chosen site as the target | technical |
| `blast_results/TEST_8_REFERENCE_blast.json` | BLAST | deletion is the owner's decision | owner | delete or keep as labeled | environmental |
| License (`pyproject.toml`: "decision pending") | project | owner's decision | owner | choose a license | environmental |
| G·U pairing, optional RNAfold (P2-5) | RNA | not implemented; RNAfold is an external tool | dependency decision (1.3.1) | implement G·U; decide on RNAfold | scientific |
| Merging the surviving primer generations (P2-1) | primer | behavior-preserving refactoring of its own | none | merge with bytecode and golden checks | technical |
| NEBcutter comparison of a circular digest (P1-4) | restriction | web service; largely covered by the Biopython pUC19 comparison | approval | decide whether still needed | environmental |
| Local BLAST on a real database | BLAST | `nt` is too large; only test databases used | approval for a ~0.5 GB download | swissprot as a proportionate step | environmental |
| HGVS conformance against a validator | mutation | no validator installed | `hgvs` package or Mutalyzer (approval) | compare a variant set | scientific |
| `core/pipeline.py` prints from `core/` | pipeline | layering follow-up of Stage F | none | move `print_status` to `ui/` | technical |
| "complexity" field names | sequence, RNA, CRISPR, mutation | kept for compatibility; meaning documented | none | rename with a documented snapshot change, if wanted | technical |
| 112 pyright errors (standard mode, Python 3.10) | primer, CRISPR, mutation, tests | reduced from 134, not cleared | none | triage for reachable `None` paths | technical |
| Tool and data ledger (1.3.1, 11.8) | project | versions recorded; source, access date, license and checksum are not | none | write the ledger | environmental |
| Line coverage | tests | no coverage tool installed | approval to install `coverage` | measure and report | technical |
| Missing-file and runtime invalid-selection tests (11.4) | input, menus | not written | none | add both | technical |
| Four references in `METHODS.md` not verified at the source | documentation | written in Stage H/I without access | approval (PubMed) | verify or keep the mark | scientific |
| `ResourceWarning` in a literature retry test | tests | a simulated `HTTPError` is not closed | none | close it in the test | technical |

---

## K. Final quality assessment

- **Scientific reliability**: good for what is implemented. Every
  calculation that has a reference implementation agrees with it
  (translation, restriction, compatibility, Tm, protein descriptors), and
  Phase 6 corrected the errors those comparisons and the audit found.
  Heuristic modules are labeled as heuristics. Nothing is validated
  against experiments, and three comparisons that need external access
  are open.
- **Functional completeness**: every advertised feature works or is
  stated as not implemented; the larger absent features (planner, CFD,
  Type IIS, genome-wide specificity, 5'-partial ORFs) are listed, not
  hidden.
- **Architecture**: much improved and consistent at the boundaries (one
  input contract, a complete API, one presentation layer). The analysis
  modules are still large files that mix menus and science.
- **Maintainability**: improved by the removal of the primer generations,
  the shared helpers, enforcement tests and documentation that matches
  the code. The size of individual modules (up to ~5,600 lines) remains a
  cost.
- **Terminal interface consistency**: high and enforced; three remaining
  deviations are documented.
- **Language quality**: consistent terminology, spelling and message
  conventions, enforced where mechanical; the remaining judgment lies in
  long scientific explanations, which were reviewed by inventory and
  sampling rather than line by line.
- **Reproducibility**: strong within the environment: hashed inputs and
  outputs, versioned reports, recorded parameters, a network-free suite.
  Weak across environments: one Python version and one operating system.
- **Testing maturity**: broad (unit, golden master, reference,
  integration, interface and style enforcement). Passing tests show that
  the implementation matches its specification and its references; they
  are not independent scientific validation.

Synapse is not finished software and is not presented as such. It is a
coherent, honest and well-tested basis whose remaining gaps are known,
listed, and in most cases a matter of external resources or the owner's
decisions rather than of undiscovered defects.
