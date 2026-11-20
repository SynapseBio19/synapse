# Synapse change log (scientific and behavioral changes)

Every change that alters a scientific result or a default is recorded here
with the reason, the old and new behavior, and the evidence. Golden-master
snapshots (`tests/golden/`) are regenerated only together with an entry here.

Verification levels used below:

- **Software correctness** – the code does what its specification says (unit tests).
- **Method conformity** – the result agrees with the published method or an
  independent reference implementation (known-answer / cross-validation tests).
- **Biological validity** – the result is confirmed by experiment or curated
  data. Synapse does not claim this for any output.

## Stage 0 – baseline (2026-10-29)

- Git repository initialized; tag `baseline-audit` = project exactly as audited.
- Added `tests/` (standard-library `unittest`, no extra dependencies):
  - `test_builtin_selftests.py` runs the 13 built-in module self-tests.
  - `test_golden_master.py` + `golden/*.json`: snapshots of current behavior
    of all 12 analysis modules.
  - `test_p0_regressions.py`: known-answer tests for the nine P0 audit
    findings, initially marked `expectedFailure` (20 tests reproduce the defects).
- Added `pyproject.toml` (Python >= 3.10, no runtime dependencies).
- No source file changed.

### Compatibility fix: `analyses/protein.py` (Python 3.10/3.11)

- Problem: seven f-strings in report-printing functions spanned several lines
  inside `{...}`. That syntax is only valid from Python 3.12 (PEP 701); on
  3.10/3.11 the whole protein module failed to import. Found by running the
  new test suite on Python 3.10.
- Change: the expressions were joined onto one line. Pure syntax change, no
  behavior change (golden master unchanged).

## Stage 1 – correction of the nine P0 audit findings (2026-10-29)

All 20 P0 known-answer tests in `tests/test_p0_regressions.py` (plus 4 new
edge-case tests) now pass. Golden-master snapshots for blast, cloning, orf,
restriction and variant were regenerated; every difference was reviewed and
is listed below. Snapshots for sequence, rna, protein, mutation, crispr,
literature and primer did not change.

### P0-1 BLAST query coverage and evidence band (`analyses/blast.py`)

- Problem: `BlastHSP.query_coverage` returned the non-gap fraction of the
  alignment, i.e. 100 % for every gap-free HSP. The stored run
  `blast_results/TEST_8_REFERENCE_blast.json` reported 100 % for a 32-bp HSP on
  a 99-bp query. Evidence bands ignored significance: E = 7.45 was labeled
  "MODERATE_SEQUENCE_EVIDENCE".
- Now: HSP coverage = HSP query interval / query length; hit coverage
  (`best_coverage_percent`) = union of all HSP query intervals / query length,
  as NCBI's "Query cover". The old value is kept as `alignment_nongap_percent`.
  Evidence bands require best E-value <= `DEFAULT_EVIDENCE_EVALUE_THRESHOLD`
  (1e-3), otherwise `NOT_SIGNIFICANT`. The threshold is a parameter
  (`significance_evalue`); 1e-3 is a common convention, not a universal value
  (E-values depend on database size and query length).
- Evidence: TEST_8 recomputes to 32.3 % coverage; golden scenario 100 % -> 20 %.
- Built-in self-test fixture `_test_evidence_profile` used a 4-nt query with a
  20-nt HSP (internally inconsistent); it now uses the matching 20-nt query.
- Not yet validated against live NCBI output (needs network, Stage 3 plan).

### P0-2 ORF hand-over to CRISPR and BLAST (`core/orf_context.py`, `analyses/crispr.py`, `analyses/blast.py`)

- Problem: both modules read the 0-based ORF `start` as 1-based and ignored the
  strand; 42/42 CRISPR guide coordinates were wrong for a selected ORF.
- Now: one shared adapter `core.orf_context.resolve_orf_region()` locates the
  ORF on the parent sequence, checks that the sequence matches, and returns
  1-based forward-strand coordinates and the strand. CRISPR designs on the
  ORF's forward-strand genomic region, so guides on both strands carry parent
  coordinates. BLAST keeps the coding-orientation query and records the
  1-based region and strand. Generic ORF-like dicts without Synapse fields keep
  the previous 1-based interpretation (built-in self-tests unchanged).
- Also fixed: nested "Select ORF / Region" in the CRISPR menu lost the region
  offset.

### P0-3 Overhang polarity in end compatibility (`analyses/cloning.py`, `analyses/restriction.py`)

- Problem: 5' and 3' overhangs with identical letters were declared
  compatible (SacI + HindIII, SphI + NcoI).
- Now: `get_overhang_polarity()`; ends must share polarity and sequence.
  Golden changes: HindIII<->SacI and NcoI<->SphI now incompatible.

### P0-4 Enzyme data (`data/restriction_enzymes.py`, `analyses/cloning.py`)

- ClaI corrected from A^TCGAT (4-nt overhang) to AT^CGAT (2-nt 5' overhang CG).
- Cloning overhangs are now derived from the database cut positions instead
  of a hand-written table that had NdeI "CATATG" (correct: TA) and ClaI "TCGA".
- Golden changes: ClaI now compatible with TaqI/MspI/HpaII (all 5'-CG) and no
  longer with XhoI/SalI (5'-TCGA). These agree with the published compatible-
  cohesive-end relationships for these enzymes.
- Added `methylation_note` for DpnI, ClaI, XbaI, HpaII (display only; digests
  still do not model methylation). `ENZYME_DATABASE_VERSION` added.

### P0-5 Ambiguity-dependent sites in digests (`analyses/restriction.py`)

- Problem: POSSIBLE sites (e.g. inside a run of N) were counted as cuts;
  10 N in 110 bp gave 17 cuts and six 1-bp fragments.
- Now: only CONFIRMED sites cut by default; POSSIBLE sites are listed in
  `possible_cuts` and in the UI note; `include_possible=True` gives the old
  worst-case view.

### P0-6 Start-loss classification (`analyses/mutation.py`)

- Problem: every Met->X change was "Start-loss".
- Now: `classify_coding_substitution(..., codon_index=)`; Start-loss only for
  codon 1. Limitation: alternative initiation codons in the reference (GTG/TTG)
  are not recognized as starts by the mutation module.

### P0-7 Variant normalization (`analyses/variant_analysis.py`)

- Problem: with a reference sequence, a non-repeat deletion `7 AG>A` was
  rewritten to `4 AA>A` (a different variant); shifting stopped one base too
  late (no anchor); REF was never checked.
- Now: Tan, Abecasis & Kang (2015) Algorithm 1 (right-trim / left-extend loop,
  then left trim); VCF right-anchoring at position 1; REF is verified and a
  mismatch adds `annotation["validation_flags"] = ["REF_MISMATCH"]` and
  disables shifting. Symbolic, `*` and non-ACGTN alleles are left unchanged.
- Not yet cross-checked against `bcftools norm` (tool not available here).

### P0-8 Initiator methionine (`analyses/orf.py`, `analyses/protein.py`, `core/pipeline.py`)

- Problem: ORFs starting with GTG/TTG were translated with V/L at position 1.
- Now: `translate_orf_sequence(..., initiator_start=True)` and
  `analyze_translation(..., initiator_start=False)`; protein analysis started
  from a selected ORF uses `initiator_start=True`. Arbitrary sequence windows
  keep the elongation code (GUG stays Val).

### P0-9 ORF defaults (`analyses/orf.py`) – configurable

- Problem: minimum length 10 nt, ATG/GTG/TTG starts for every organism, every
  nested start reported and a heuristic ranking put a 36-bp TTG ORF above the
  99-bp ATG CDS of TEST_8_REFERENCE.
- New defaults, following NCBI ORFfinder conventions: minimum 75 nt (stop
  included), ATG only, longest ORF per stop codon, ranking by length.
  These are screening conventions, not biological thresholds.
- Alternatives remain available: parameters `start_codon_policy`,
  `nested`, `ranking`, `min_length`, the new ORF menu "[6] Search Settings",
  and `LEGACY_ORF_SETTINGS`, which reproduces the previous output exactly
  (test `test_legacy_settings_reproduce_baseline_snapshot`, snapshot
  `tests/golden/orf_legacy_baseline.json`).
- Search parameters are stored in the analysis result and in the selected ORF
  (`search_parameters`). The score formerly shown as "Confidence" is shown as
  "Heuristic score" (it is not a probability).
- Visible consequence: short ORFs (< 75 nt, e.g. TEST_8_SNP_NONSENSE, 45 nt)
  are no longer listed unless the minimum length is lowered.

## Stage 2 – shared data model, visible result storage, input, logging (2026-10-30)

No scientific algorithm changed in this stage; all golden snapshots are
unchanged.

- `core/models.py`: `Interval` (0-based half-open + strand, with explicit
  1-based conversions), `Feature`, `SequenceRecord` (sequence + provenance:
  file path, format, file SHA-256, sequence SHA-256, topology) and
  `AnalysisResult` (module, data, parameters, input hash, warnings,
  limitations, timestamp, Synapse and Python version).
- `core/seqtools.py`: shared reverse complement / IUPAC helpers and NCBI
  genetic codes 1, 2, 4 and 11 in NCBI's compact table format. Table 1 is
  verified against all four independently written codon tables already in
  the project (test). Not yet used by the older modules (Stage 4).
- `core/pipeline.py`: generic result store (`store`, `store_result`,
  `latest`, `results`), `store_primer_analysis`, provenance via
  `attach_record`; changing the selected ORF invalidates ORF-dependent
  results (RNA, protein, CRISPR, BLAST). Status view lists all 12 modules,
  the sequence SHA-256 and the source file.
- Results that were previously discarded are now stored: primer, mutation,
  variant (wrong method/keys before), protein (no pipeline passed), CRISPR
  and BLAST (menus now return their last result to `main.py`).
- Silent `except: pass` around result storage replaced by logged warnings.
- `fileio/loader.py` + `fileio/genbank.py`: load FASTA or GenBank chosen by
  the user (previously only `examples/dna.fasta`). GenBank topology from
  LOCUS; FASTA `topology=circular` header token; file SHA-256 recorded.
  GenBank location parsing covers simple spans, complement and join; other
  forms keep only the raw string (documented limitation).
- `core/logging_setup.py`: rotating log file `logs/synapse.log`; `main.py`
  errors show a short message and write the traceback to the log.
- UI: sequence list shows length and topology; `[F]` loads another file;
  `[P]` shows pipeline status; ORF confirmation shows the 1-based location
  and strand; duplicated "RNA CONNECTED" message removed from `main.py`.
- Tests: `tests/test_stage2_foundation.py` (15 tests) and an end-to-end
  menu test `tests/test_main_smoke.py` (ORF -> RNA -> protein -> status).

## Stage 3 – P1 scientific improvements and reproducibility (2026-10-30)

### Primer melting temperature (`analyses/primer.py`) – output changes

- Default nearest-neighbor parameters: SantaLucia (1998) unified set with
  the terminal-AT/GC initiation terms, the 0.368·(N−1)·ln[Na⁺eq] entropy salt
  correction, C/4 (C for self-complementary oligos) and the −1.4 e.u.
  symmetry term. This is the model documented for Primer3
  (`PRIMER_TM_FORMULA=1`, `PRIMER_SALT_CORRECTIONS=1`).
- Divalent cations: [Na⁺eq] = [Na⁺] + 120·√([Mg²⁺] − [dNTP]) (von Ahsen et al.
  2001), replacing the previous ad-hoc rule [Na⁺] + 4·([Mg²⁺] − [dNTP]).
- Selectable alternatives: `nn_parameter_set` = `santalucia_1998` (default),
  `santalucia_2004`, `synapse_legacy`; `salt_model` =
  `monovalent_equivalent` (default, von Ahsen), `na_only`,
  `synapse_legacy_4x`. `synapse_legacy` + `synapse_legacy_4x` reproduce the
  previous values to 5 decimals (test).
- Verification: Primer3 manual worked example (CGTGACGTGACGGACT, 50 mM Na⁺,
  50 nM): ΔH −128.8 kcal/mol and salt-corrected ΔS −361.736 cal/(mol·K)
  reproduced exactly; Tm 50.550 °C vs. 50.554 °C in the manual (the manual
  uses R = 1.987). An independent re-implementation of the published
  equations in the test file agrees to 1e-6 °C for 6 primers × 4 buffer
  conditions. Running the Primer3 or Biopython executables was not possible
  in this environment, so no cross-tool run is claimed.
- Visible consequence (golden `primer.json` regenerated): with 1.5 mM Mg²⁺ /
  0.2 mM dNTP the calculated Tm is ~6–7 °C higher than before, so the
  designer selects shorter primers (18–20 nt instead of 20–24 nt) for the
  same 60 °C target. Neither value is an experimental measurement; real
  annealing temperatures must be determined empirically.

### Variant analysis (`analyses/variant_analysis.py`)

- Multi-allelic splitting reads `Number=` from the `##INFO`/`##FORMAT`
  header (VCF 4.3 reserved keys as fallback): `A` keeps the selected ALT,
  `R` keeps [REF, ALT], `G` (PL/GL/GP) keeps genotypes 0/0, 0/a, a/a using
  F(j,k) = k(k+1)/2 + j. Previously only list length was used.
- GT of other ALT alleles is set to 0 (convention of `bcftools norm -m-`).
- `.vcf.gz` (gzip/BGZF) files are read; the UI accepts a VCF file path.
- FILTER "." now counts as `unfiltered` and no longer as `failed_filters`
  (VCF: "." = filters not applied). Golden `variant.json` regenerated; the
  only differences are `failed_filters` 1 → 0 and the new `unfiltered: 1`.
- `check_reference_alleles()` compares REF with a supplied reference per
  contig and flags `REF_MISMATCH`; `normalize_variants(contig_references=)`
  normalizes each contig against its own sequence.
- Not done: a comparison against `bcftools norm` output (tool not available
  here). Recommended validation for a later phase.

### Restriction analysis (`analyses/restriction.py`)

- `topology="circular"` finds sites spanning the origin (`wraps_origin`)
  and maps cut positions modulo the sequence length. Linear behavior is
  unchanged (golden unchanged).

### BLAST (`analyses/blast.py`)

- Search parameters reported by BLAST in the XML output are stored
  (`reported_parameters`), plus database provenance (`database_details`;
  local databases: title, date, version, counts from `blastdbcmd -info`).
- Remote searches ask for confirmation before the query sequence is sent to
  NCBI (privacy). Nothing changes for local searches.
- Not done: resuming a remote search from a saved RID (deferred).
- Local BLAST against `data/blast_db/` could not be run (databases are
  intentionally excluded from this copy); only parsing is tested.

### Literature (`analyses/literature.py`)

- Per-host request throttling (NCBI: 3 req/s, 10 req/s with
  `SYNAPSE_NCBI_API_KEY`) and up to 3 attempts on HTTP 429/5xx, honoring
  `Retry-After`.
- Partial dates are expanded correctly (YYYY-MM → last day of the month for
  upper bounds). Bug fixed: OpenAlex received `YYYY-MM-01-01`.
- PubMed date filters use `YYYY/MM/DD`; optional `api_key`.
- When opened from the main menu, the sequence header is offered as a
  search suggestion (shown and only used after the user presses Enter); no
  sequence data are sent.

### ORF genetic code (`analyses/orf.py`) – P1-9

- New option `genetic_code` (NCBI tables 1, 2, 4, 11; menu [8] in ORF search
  settings). Default 1, so existing results are unchanged (golden
  unchanged). Stop codons, translation, ambiguity and codon usage follow the
  selected table; `atg_only` keeps ATG, `atg_and_alternative` uses the
  NCBI start codons of the selected table; new policy `table_starts`.
- ORFs carry `genetic_code`; parameters record `genetic_code_table`.
- Limitation: the RNA and protein modules still translate with table 1;
  `main.py` prints a warning when an ORF from another table is selected.
  5′-partial ORFs and ORFs across the origin of circular sequences are not
  implemented (deferred).

### Mutation analysis (`analyses/mutation.py`) – P1-10

- The coding-consequence step can take the CDS from the selected pipeline
  ORF (both strands; the variant CDS is derived only when the variant has
  the reference length).
- `assess_coding_frame()` warns when the reference CDS is not a multiple of
  3, does not start with a start codon, contains in-frame stops or lacks a
  terminal stop; the frame assumption is recorded in the result notes.
- Comparisons longer than 2,000 nt get a note (O(n·m) global alignment).
- `mutation_hgvs_c/p` are documented as HGVS-like and not validated
  against the HGVS specification.

### Reproducible analysis report (`core/report.py`, menu `[R]`) – P1-12

- Exports all results stored in the session as JSON (complete payloads) or
  Markdown (summary): input provenance and sequence SHA-256, per analysis
  the parameters, warnings, limitations, timestamp, versions and a SHA-256
  of the canonical JSON payload.

### Tests

- `tests/test_primer_thermodynamics.py` (7), `tests/test_stage3_additions.py`
  (21) and a menu test for `[R]` in `tests/test_main_smoke.py`.

## Stage 4 – consolidation and preparation for a future UI (2026-10-30)

No scientific result changes in this stage; all golden snapshots unchanged.

- **Primer patch chain**: 13 superseded function definitions in
  `analyses/primer.py` (1,014 lines; earlier versions of e.g.
  `design_primers_professional`, `generate_professional_primer_candidates`,
  `run_primer_analysis_ui`) were removed. They were never reachable: a later
  definition replaced the name before any call. Evidence of equivalence:
  the compiled bytecode of every remaining module-level function and class
  is identical before and after (automated comparison), and all tests pass.
  The live v5–v11 wrapper chain for candidate ranking was kept; flattening
  it would change code that is actually executed and is left for later.
- **Shared tables**: the standard genetic code in `orf`, `mutation`, `rna`
  and `data/genetic_code.py` and the IUPAC complement tables in `crispr`,
  `primer` and `mutation` are now derived from `core/seqtools.py`
  (`codon_dict`). Each replaced table was asserted equal to the derived one
  before replacement. One former hand-written table is kept as an
  independent test fixture (`tests/fixtures/standard_code_reference.py`).
  The modules' own normalization and error handling were not changed.
- **Genetic code in protein analysis**: `analyze_translation`,
  `analyze_complete_protein`, the protein menu and the pipeline use the
  genetic code of the selected ORF (default table 1).
- **`core/api.py`**: non-interactive entry points (load, sequence overview,
  ORF, RNA, protein, restriction digest, primer design, sequence comparison,
  VCF, CRISPR) returning `AnalysisResult`. The ORF entry point is tested to
  give results identical to the module call.
- **Documentation**: `docs/METHODS.md` (methods, defaults, references,
  validation status per module); README extended.
- Tests: `tests/test_stage4_consolidation.py` (10).

## Stage 5 – unresolved residues keep their position (2026-10-30)

- `analyses/protein.py`: new `get_positional_protein()` keeps unresolved
  residues ("?") instead of deleting them. The position-dependent analyses
  (hydropathy profile, hydrophobic regions, transmembrane screen,
  flexibility profile, signal-peptide screen, instability index) use it and
  skip only the windows or dipeptides that contain a "?". Previously such
  residues were removed, which silently shifted every later position: for a
  75-residue protein with one ambiguous codon, every reported window
  position after that residue was off by one.
- Visible consequence only for sequences that contain ambiguous codons:
  windows covering an unresolved residue are no longer reported (instead of
  being reported with a wrong position), and the instability index is
  normalized by the number of scored dipeptides. Fully resolved sequences
  give identical values; all golden snapshots are unchanged.
- Tests: `tests/test_stage5_protein_positions.py` (5), including the
  position-identity check against the fully resolved sequence.

## Assessment against the audit (2026-10-30)

`docs/AUDIT_STATUS.md` lists every P0, P1, P2 and P3 recommendation with one
of: implemented and tested / implemented but not independently validated /
partial / deferred / rejected, including the validation steps that could not
be performed in this environment.

## Phase 6 – Tier 1 corrections (2026-10-31)

Defect IDs refer to `docs/PHASE6_INITIAL_AUDIT.md`. Each defect is corrected
in its own commit with its own regression tests in
`tests/test_phase6_tier1.py`.

### D12 Contig names in REF checking and normalization (`analyses/variant_analysis.py`, `core/api.py`)

- Problem: reference sequences were looked up by exact CHROM name, and a
  name that did not match was skipped silently. With a VCF using `1` and a
  reference keyed `chr1` (or the reverse), `check_reference_alleles`
  reported `mismatched: 0` (the skip was visible only in `not_checked`),
  `normalize_variants` returned the variant trimmed but not left-aligned
  and without any flag, and `core.api.analyze_vcf` recorded
  `parameters["normalized"] = True`. Example: `1 7 AA>A` against
  `chr1 = GGGAAAAAGGG` stayed `7 AA>A`; the correct left-aligned form is
  `3 GA>G`.
- Now:
  - `resolve_contig_names()` maps CHROM names onto reference names. Mode
    `contig_matching="exact"` (default) keeps exact matching. Mode
    `"chr_prefix"` (opt-in) also pairs a nuclear contig with the name that
    differs only by a leading `chr`; an exact name always wins.
  - Mitochondrial names (`M`, `MT`, `chrM`, `chrMT`) are never paired:
    hg19 chrM (16,571 bp) and the rCRS MT of GRCh37/GRCh38 (16,569 bp) are
    different sequences. Case variants (`Chr1`) are never folded. CHROM is
    never renamed.
  - A variant without a reference is flagged `REFERENCE_CONTIG_NOT_FOUND`
    with a note that it was not REF-checked or left-aligned; in exact mode
    the note names the `chr`-prefix partner when one exists. A variant
    matched through the prefix rule is flagged `CONTIG_NAME_ALIASED` and
    records `annotation["reference_contig"]`.
  - `check_reference_alleles()` additionally reports `contigs_not_found`,
    `aliased_contigs`, `suggested_aliases`, `contig_matching`, and splits
    `not_checked` into `not_checked_no_reference` and
    `not_checked_symbolic_or_non_acgtn`. Existing keys are unchanged.
  - `core.api.analyze_vcf(contig_matching=)`: `parameters["normalized"]` is
    True only when every variant had a reference; `data["normalization"]`
    lists the contigs without one and the counts; each missing contig and
    each alias produces a warning.
  - Variant menu `[8] Normalize Variants` reports how many variants had a
    reference sequence and how many were only trimmed, instead of
    "Normalized N variant(s)".
  - `normalize_variants()` no longer writes its new flags into the input
    variants (copy-on-write of `validation_flags`).
- Design source: guidance of the project-local `genomic-coordinates` skill
  (`references/reference-builds.md`, "Contig naming"): renaming between
  naming schemes is valid only when the sequences are identical, and
  `chrM`/`MT` must not be mapped silently. The prefix rule is therefore
  opt-in, recorded per variant, and still subject to the per-variant REF
  check against the supplied sequence.
- Unchanged: exact-name behavior, left-alignment, trimming, every golden
  snapshot (no snapshot regenerated), the variant self-test (27/27).
- Tests: 9 new tests (13 cases with sub-tests). Suite 109 → 118, all
  passing. pyright at `pythonVersion = 3.10` reports the same 5
  pre-existing diagnostics in the changed files before and after; none is
  on a changed line.
- Not done: `bcftools norm` comparison (no native Windows build); contig
  length checks against `##contig` header lines; RefSeq/GenBank accession
  naming (`NC_000001.11`), which needs an assembly report rather than a
  string rule.

### D4 Cloning used its own site finder without the CONFIRMED/POSSIBLE distinction (`analyses/cloning.py`, `analyses/restriction.py`)

- Problem: `cloning.py` had a private recognition matcher
  (`_recognition_matches`) that returned every IUPAC-compatible window as a
  site. The P0-5 fix (ambiguity-dependent sites do not cut) existed only in
  `restriction.py`. `cloning.find_restriction_sites("A"*10 + "N"*6 + "A"*10,
  "GAATTC")` returned `[11]`; a vector whose only EcoRI-compatible site was
  `GANTTC` gave `evaluate_cloning_strategy(...) == "VALID"`, and the
  automatic planner reported a valid strategy built on that site.
- Now:
  - One matcher: `restriction.find_recognition_site_matches(sequence,
    recognition_sequence, topology)` holds the existing loop of
    `find_restriction_site_matches`, which now delegates to it (no behavior
    change; the pUC19 cross-check against Biopython still gives 0
    differences in 132 comparisons). `cloning.find_restriction_sites` and
    `cloning.analyze_enzyme` use it; cloning's private matcher was removed
    and `cloning.bases_match` is the restriction function.
  - `cloning.analyze_enzyme`: `sites` / `site_count` hold CONFIRMED sites
    only; new `possible_sites`, `possible_site_count`, `site_matches`.
    `find_restriction_sites(..., include_possible=False)`.
  - An enzyme with a POSSIBLE site in a molecule is neither a confirmed
    unique cutter nor a confirmed non-cutter there — the rule
    `restriction.py` already applies in `find_unique_cutters`,
    `recommend_enzyme_pairs` and `get_cloning_compatible_sites`. Cloning
    ends (`_end_definition`), single-enzyme analysis, enzyme suitability and
    the planner diagnostics (`AMBIGUITY_DEPENDENT`) apply it and state the
    reason with the positions.
- Interface: signatures are unchanged except for the new
  `include_possible` keyword; `cloning.find_restriction_sites` and
  `analyze_enzyme` now raise `ValueError` for an empty or non-IUPAC
  sequence instead of silently skipping invalid windows. The planner and
  both strategy evaluators validate their input first and are unaffected.
- Unchanged: results for sequences without ambiguity symbols (golden
  snapshots unchanged; cloning 13/13, cloning UI and restriction 14/14
  self-tests pass). Planner on a 1.3 kb vector, 10 enzymes, same process:
  old 100–110 s, new 68–81 s, identical results.
- Tests: 7 new (`CloningAmbiguousSiteTests`), including identity of the
  cloning and restriction site lists for all 22 enzymes on a sequence with
  IUPAC symbols. Suite 118 → 125.
- Found and recorded, not changed: S17 (the cloning digest uses the
  recognition start instead of the cut coordinate) and S18 (planner
  runtime); see `docs/PHASE6_INITIAL_AUDIT.md` §5.

### D3 Restriction map and site-based results ignored circular topology (`analyses/restriction.py`)

- Problem: only `perform_digest` honored `topology="circular"`. The
  restriction map, per-enzyme results, unique cutters, non-cutters,
  internal-site warnings and enzyme-pair recommendations searched linearly,
  so `analyze_restriction(..., "circular")` on a plasmid whose EcoRI site
  spans the origin returned a digest that cuts once together with a map
  without the site and EcoRI listed as a non-cutter. The map menu offered
  no topology choice at all.
- Now: `topology` (default `"linear"`) is accepted and passed through by
  `find_restriction_sites`, `analyze_restriction_site`, `analyze_enzyme`,
  `create_restriction_map`, `summarize_restriction_map`,
  `find_unique_cutters`, `find_non_cutters`, `find_internal_sites`,
  `recommend_enzyme_pairs`, `summarize_enzyme_sites` and the `get_*`
  wrappers; `analyze_restriction` passes it everywhere.
  `get_cloning_compatible_sites` takes `vector_topology` and
  `insert_topology`. For an origin-spanning site, end and cut coordinates
  are reported modulo the sequence length (the digest's convention) and
  the overhang is read across the origin; map entries carry `wraps_origin`.
  The Restriction Map menu asks for the topology like the other restriction
  menus, defaulting to the loaded sequence's topology.
- Unchanged: every linear result, and the keys of `analyze_restriction_site`
  (part of the golden master); golden snapshots unchanged.
- Validation: the circular map agrees with Biopython 1.88 `Bio.Restriction`
  on pUC19 (L09137.2) in 132/132 comparisons (22 enzymes × 2 topologies ×
  3 rotations, including a site across the origin and a cut exactly at
  the origin).
- Tests: 7 new (`RestrictionCircularTopologyTests`), including a check
  that map and digest cut sets and the non-cutter list agree for all 22
  enzymes in both topologies with a site rotated across the origin. Suite
  125 → 132.

### D5 Cloning analyzed circular vectors as linear (`analyses/cloning.py`, `core/pipeline.py`, `analyses/restriction.py`)

- Problem: no cloning function accepted a topology; vector and insert
  digests were hard-coded `"linear"`. A recognition site spanning the
  origin of a circular plasmid was invisible, so a vector whose EcoRI site
  spans the origin gave `INVALID_VECTOR_ENDS`, also through
  `SynapsePipeline.run_cloning` on a sequence loaded as circular.
- Now: `vector_topology` and `insert_topology` (default `"linear"`;
  `ValueError` for other values) are accepted and passed down to site
  finding by `evaluate_cloning_strategy`, `analyze_enzyme_pair`,
  `analyze_cloning`, `assess_enzyme_suitability`, `find_compatible_enzymes`,
  `find_best_cloning_strategies`, `automatic_cloning_planner`,
  `run_advanced_cloning_analysis`, `run_cloning_analysis`,
  `get_best_cloning_strategy`, and internally by the end definitions,
  internal-site checks, digests and planner diagnostics.
  `find_restriction_sites`, `analyze_enzyme` and `build_restriction_map`
  take `topology`. Results record the topology (`vector.topology`,
  `insert.topology`, planner `topology`); vector and insert digests use it.
- Scientific convention, stated in the result (`assumptions`): on a
  circular molecule the segment between the two selected sites in
  coordinate order is taken as the region replaced (vector) or excised
  (insert). Which arc is intended cannot be inferred from the sequence.
- Integration: `SynapsePipeline.run_cloning` / `run_cloning_planner` pass
  the loaded record's topology for the vector; the Restriction & Cloning
  menu passes the sequence's topology to the cloning menu, and an insert
  chosen from the loaded records keeps that record's topology. The cloning
  menu shows both topologies.
- Unchanged: all linear results (golden snapshots unchanged; cloning,
  cloning-UI and restriction self-tests pass).
- Tests: 7 new (`CloningCircularTopologyTests`), including a circular
  insert donor, the pipeline path, and a scripted run of the menu from
  Restriction & Cloning to the planner. Suite 132 → 139.
- Not changed (recorded as S17): the cloning digest's cut coordinates.

### D1 One initiation-codon rule, taken from the genetic code (`analyses/orf.py`, `analyses/protein.py`, `core/api.py`)

- Problem: the two translators decided differently whether the first
  codon is decoded as initiator Met, and neither followed the genetic
  code. `orf.translate_orf_sequence` used a hard-coded set ATG/GTG/TTG for
  table 1, so a CTG-initiated ORF (found with `table_starts`) translated as
  `LAAA…`. `protein.analyze_translation(initiator_start=True)` — the path
  of the ORF → Protein menu and of `core.api.analyze_protein_from_cds` —
  turned any resolved first codon into Met, so a CDS beginning with GCT
  (Ala) translated as `MKF…`. The same CTG ORF therefore gave `L…` in one
  path and `M…` in the other.
- Now: `orf.initiation_codons(table)` is the union of the start codons of
  every Synapse start-codon policy for the table — the NCBI start codons
  (`core.seqtools`, verified against Biopython) plus, for table 1 only,
  the documented `atg_and_alternative` policy's GTG. Table 1:
  ATG/CTG/GTG/TTG; tables 2, 4, 11: exactly NCBI's start set.
  `orf.is_initiation_codon(codon, table)` is true when every IUPAC
  reading is an initiation codon — the rule the ORF finder already
  applies to a "definite" start. Both translators use it for the first
  codon; every other codon is read with the elongation code of the table.
  `core.api.analyze_protein_from_cds` warns when `initiator_start=True`
  but the first codon is not an initiation codon.
- Behavior changes: CTG (table 1) and definite IUPAC initiators such as
  WTG/YTG/NTG now give Met in the ORF translator; a non-initiation first
  codon is no longer forced to Met by the protein module. Unchanged: every
  ORF found by the ORF finder (its first codon is always an initiation
  codon), GTG/TTG → Met in table 1 (P0-8), all internal codons, the golden
  snapshots.
- GTG in table 1 is not an NCBI table-1 start codon; it stays an
  initiation codon because Synapse's user-selectable `atg_and_alternative`
  policy (approved with P0-8) searches for GTG-initiated ORFs. This is
  the one intentional difference from Biopython's `cds=True` rule and is
  asserted as such in the test.
- Validation: against Biopython 1.88 for tables 1, 2, 4, 11 — every
  internal codon translates as in `Bio.Seq.translate`, and the set of
  first codons decoded as Met equals Biopython's start codons (plus GTG
  for table 1) together with the codons that encode Met anyway.
- Tests: 7 new (`OrfInitiationCodonTests`; the Biopython comparison is
  skipped when the `reference` extra is not installed). Suite 139 → 146.
- Found and recorded, not changed: S19 (start-codon *classification* and
  the CDS-quality / heuristic scores still use the fixed ATG/GTG/TTG set,
  so a CTG ORF is classified "None" and scores lower).

### D2 An ORF's genetic code travels with the ORF into translation (`analyses/orf.py`, `core/api.py`)

- Problem: ORF dicts carry `genetic_code`, but both translators took only
  a DNA string with a silent table-1 default. A table-2 ORF
  (`ATG` + `TGA`×25 + `TAG`) translated with the default gave `"M"`
  instead of `"M" + "W"×25`, with no warning; the same applied to
  `core.api.analyze_protein_from_cds(orf["sequence"])` for ORFs returned
  by `core.api.find_orfs`.
- Now: `orf.translate_orf_sequence` and `core.api.analyze_protein_from_cds`
  accept the ORF dict itself and translate it with the ORF's own
  `genetic_code`; an explicit, different `genetic_code` raises
  `ValueError` instead of being used. `analyze_selected_orf` passes the
  ORF. A DNA string still defaults to table 1 (a string carries no table),
  so every existing call behaves as before; `genetic_code` now defaults to
  `None` meaning "the ORF's table, or 1 for a string".
- Editorial: the rewritten `analyze_protein_from_cds` docstring uses the
  project spelling "analyze".
- Tests: 5 new (`OrfGeneticCodeHandoffTests`), including agreement of the
  ORF translator and the API for every ORF found in a 3 kb sequence with
  tables 1, 2, 4 and 11 (several hundred ORFs). Suite 146 → 151.
- Not changed: the RNA module's three-frame codon view still uses table 1
  for a selected non-standard-table ORF. This is the documented P1-9
  limitation (`main.py` warns when such an ORF is selected), not part of
  D2; it is the remaining gap against "translation uses the selected ORF's
  genetic code" and is listed as open in `docs/FEATURE_STATUS.md`.

## Phase 6 – Sequence types (D14 with D13, 2026-10-31)

D14 was investigated in twelve sub-findings (F1–F12, listed under D14 in
`docs/PHASE6_INITIAL_AUDIT.md`); F12 is D13, and both are corrected
together. Binding decisions for this work: no percentage threshold is
presented as an established cut-off; a sequence whose type cannot be told
from its content needs an explicit type before type-dependent analysis;
protein FASTA records get a route to protein analysis. Tests are in
`tests/test_phase6_sequence_types.py`.

### D14-F2 One query-type vocabulary in BLAST (`analyses/blast.py`)

- Problem: the program table declared `query_type: "protein"` while every
  `BlastQuery` carries `sequence_type: "PROTEIN"`. Local blastp and
  tblastn therefore always failed with "blastp expects a protein query, but
  the query is PROTEIN.", and the BLAST menu offered no program at all for
  a protein query. Reproduced on the previous commit with a real
  `makeblastdb`/`blastp` run.
- Now: the table uses `"PROTEIN"` (`QUERY_TYPES = ("DNA", "PROTEIN")`),
  the menu lists programs with `programs_for_query_type`, and both the
  local engine and the NCBI client check the query against the program
  with `check_program_query_type`. The NCBI client now refuses a
  mismatched query *before* submitting it instead of sending it to NCBI.
- Validation: external tool execution — a 40-residue query against a
  two-record protein database built with `makeblastdb`, searched with
  BLAST+ 2.17.0 `blastp`, returns the correct target (skipped when BLAST+
  is not installed).
- Tests: 5 new (`BlastQueryTypeVocabularyTests`). Suite 151 → 156.
- Not changed here: how a query's type is inferred (F1, next package).

### Central sequence alphabets and type classifier (`core/alphabets.py`, new)

- Problem (F12/D13): 15 validators with 7 alphabet contracts, and no
  module able to say whether a sequence is DNA, RNA or protein.
- Now: one module with the IUPAC nucleotide and protein alphabets, a
  `normalize_sequence` that never converts or drops symbols, and
  `classify_sequence(sequence, declared=None)` returning the type
  (DNA/RNA/PROTEIN), a confidence (CONFIRMED/PROBABLE/AMBIGUOUS/INVALID),
  the consistent candidate types, reasons and invalid-symbol positions.
  `require_sequence_type` raises `SequenceTypeError` (a `ValueError`) for
  INVALID, AMBIGUOUS or a type the caller does not accept. The rules are
  in `docs/METHODS.md`; there is no frequency threshold.
- No analysis uses it yet; it is wired in by the following packages, so
  this commit changes no existing behavior.
- Validation: the letter sets equal Biopython's `Bio.Data.IUPACData`
  (external reference comparison, skipped without Biopython); alphabet
  sources checked in PubMed (PMIDs 2582368, 6692818).
- Tests: 11 new (`SequenceClassifierTests`). Suite 156 → 167.

### D14-F3, F11 Declared molecule type is read, kept and reported (`fileio/`, `core/models.py`, `core/report.py`)

- Problem: a GenPept record (`44 aa`) loaded as DNA, a GenBank RNA record
  failed the DNA gate, FASTA `[moltype=...]` was ignored, and
  `SequenceRecord.alphabet` was always `"DNA"`; the type appeared nowhere
  in the pipeline or the exported report.
- Now: the GenBank reader records the LOCUS molecule (`protein` for `aa`
  records, otherwise the token after `bp`, e.g. `mRNA`); the FASTA loader
  reads the NCBI `[moltype=...]` modifier (`[moltype=protein]` is a
  Synapse extension). Every loaded record is classified with
  `core.alphabets` (`alphabet`, `declared_type`, `sequence_classification`);
  a declared type the letters contradict adds a note. NCBI's own default
  of "genomic DNA" for a FASTA record without the modifier is not
  adopted: that default is exactly what let a peptide pass as DNA.
- `SequenceRecord` classifies itself: `alphabet` is DNA/RNA/PROTEIN or
  None when undetermined, with `type_confidence`, `type_reasons`,
  `declared_type`, `molecule_type`, and `with_declared_type()` for a
  user's declaration. An `alphabet` passed to the constructor counts as a
  declaration. `provenance()` and the report's input section (JSON and
  Markdown) include the sequence type and its confidence.
- Compatibility: `alphabet` was never read before; existing records of
  plain DNA classify as DNA. The main menu does not use the type yet.
- Tests: 6 new (`DeclaredMoleculeTypeTests`). Suite 167 → 173.

### D14-F1, S20 BLAST takes the query type from `core.alphabets` (`analyses/blast.py`)

- Problem: `infer_sequence_type` called a sequence DNA when every letter
  was an IUPAC nucleotide code and protein otherwise. A peptide such as
  `MKVHDANWYSK...` was searched with blastn (0 hits, no warning); RNA, a
  T/U mixture and gapped DNA were reported as invalid *protein*. The
  protein alphabet lacked J, U and O although BLAST+ accepts them (S20,
  confirmed with BLAST+ 2.17.0 blastp: U and J are accepted, O is scored
  as X with a warning).
- Now: `infer_sequence_type(sequence, declared=None)` delegates to
  `classify_sequence` and refuses INVALID and AMBIGUOUS input instead of
  guessing. An RNA query is searched as DNA with U read as T, recorded in
  the new `BlastQuery.notes`. FASTA queries honor `[moltype=...]`; an
  undetermined record without it is refused with a hint. The menu asks
  "DNA or Protein" for an undetermined selected sequence, and
  `blast_menu` takes an optional `sequence_type` so a declared type is
  not asked again. The protein alphabet is `core.alphabets.PROTEIN_IUPAC`.
- Validation: external tool execution — the IUPAC-letter peptide,
  declared a protein, finds itself with local blastp.
- Tests: 7 new (`BlastQueryTypeInferenceTests`). Suite 173 → 180.

### D13, D14-F4/F5/F7 One input contract for the menu gate and the API; protein records get their own route (`main.py`, `core/api.py`, `core/validation.py`, `core/seqtools.py`, `analyses/protein.py`, `core/pipeline.py`, `ui/menu.py`)

- Problem: two `normalize_dna` functions with different U handling (D13);
  the menu gate accepted a peptide made of IUPAC-nucleotide letters and
  the empty string (F4); `core.api` and primer design silently converted
  RNA to DNA while the menu rejected it (F5); the API checked no types,
  so `sequence_overview(protein)` reported a GC content and
  `analyze_protein_from_cds(protein)` returned `"????????"` (F7); a
  protein record could only reach the DNA overview.
- One contract: `core.validation.normalize_dna`, `core.seqtools.normalize_dna`
  and `core.alphabets.normalize_sequence` all strip whitespace and
  uppercase, and convert nothing. `validate_dna("")` is now False.
  `GeneticCode.translate` still reads U as T codon by codon.
- Menu gate (`main.py`): each record's type comes from its
  `SequenceRecord` classification. An undetermined record (e.g. IUPAC
  ambiguity codes other than N) is not guessed: the user declares DNA, RNA
  or protein, and a declaration the letters contradict is refused. An RNA
  record is analyzed in the DNA alphabet with a printed note, and the
  report records `analyzed_alphabet`. Invalid records list the offending
  symbols. A protein record opens a protein menu (protein analysis,
  BLAST with a protein query, report export) instead of the DNA overview.
  The record list shows `aa` for protein records.
- API (`core.api`): every entry point checks its input type and takes an
  optional `sequence_type`. Nucleotide analyses accept DNA and RNA; RNA is
  converted to the DNA alphabet with a warning in the result. Restriction
  digests and CRISPR guide design accept DNA only. New
  `analyze_protein(record_or_sequence)` for protein sequences. A record's
  declared type is used. Wrong or undetermined types raise
  `SequenceTypeError` (a `ValueError`).
- `SynapsePipeline.attach_record` accepts an RNA record for a sequence
  loaded in the DNA alphabet; BLAST and report export in the main menu use
  shared helpers (`store_blast_run`, `export_report_interactive`) so the
  DNA and protein menus behave the same.
- Behavior change, intended: DNA records with IUPAC ambiguity codes
  other than N now ask for the type once (e.g. `TEST_8_AMBIGUOUS` in
  `examples/dna.fasta`); a header modifier `[moltype=genomic DNA]`
  declares it in the file. Plain DNA records behave as before; the golden
  master is unchanged.
- Tests: 9 new (`InputContractTests`, `MainMenuSequenceTypeTests`, the
  latter driving `main.py` through a protein, an undetermined, an RNA and
  an invalid record). Suite 180 → 189.
- Not changed here: protein residues X/B/Z/J/U/O in the protein module
  (F6), lower-case ORF search (F8), `transcribe_rna` (F9), RNA GC content
  in the sequence module (F10).

### D14-F6 IUPAC protein codes keep their position and meaning (`analyses/protein.py`)

- Problem: X, B, Z, J, U (selenocysteine) and O (pyrrolysine) were
  "invalid"; the positional view deleted them, so every later sliding
  window and dipeptide shifted (the Stage 5 defect, live for user-supplied
  proteins); U and O had no mass; an internal `*` was deleted the same
  way.
- Now: every position is kept. In position-dependent analyses every
  non-standard symbol becomes `?`, and the windows containing it are
  skipped, as for unresolved translation positions; one terminal `*` is a
  stop and is removed. Quality: B, Z, J, X count as ambiguous, U and O as
  `nonstandard_residues` (valid), an internal `*` as invalid, and
  `terminal_stop` is reported. Composition reports `nonstandard_counts`.
  Molecular weight: U and O use residue masses derived from Biopython's
  `IUPACData.protein_weights`; B, Z, J give a range over their residues
  (J: I and L have the same mass); X and `?` leave it undetermined, as
  before. The quality printouts show a "Sec/Pyl (U, O)" line only when
  such residues occur.
- Test changed: `tests/test_stage5_protein_positions.py` asserted
  `get_positional_protein("MK?L*x") == "MK?L"`, i.e. the deletion. It now
  expects `"MK?L??"` (positions kept) and `"MK?L*" -> "MK?L"` (terminal
  stop). The new expectation is stricter, not weaker.
- Limitation: hydropathy, charge and pI scales have no values for U and O;
  they are treated as unresolved there.
- Validation: U/O residue masses equal Biopython's within 0.001 Da, and
  molecular weights of U/O-containing peptides agree with
  `Bio.SeqUtils.molecular_weight` within 0.1 Da (the standard-residue
  tables differ slightly; external reference comparison, skipped without
  Biopython).
- Tests: 5 new (`ProteinResidueCodeTests`). Suite 189 → 194. Golden master
  unchanged.

### D14-F8, F9, F10 Module entry points normalize and check their input (`analyses/orf.py`, `analyses/rna.py`, `analyses/sequence.py`)

- F8: `orf.find_orfs` on lower-case (soft-masked) DNA found 0 ORFs, and
  RNA or protein input also gave a silent 0. It now uppercases and strips
  whitespace, and raises `ValueError` naming the first offending symbols
  for anything outside the IUPAC DNA alphabet (RNA included; `core.api`
  converts RNA explicitly).
- F9: `rna.transcribe_rna` transcribed anything (`MKTAY` → `MKUAY`). It
  now normalizes, returns RNA unchanged, and raises for non-nucleotide
  symbols.
- F10: `sequence.count_nucleotides` ignored U, so the GC content of an RNA
  passed directly to the DNA module omitted U from the denominator
  (64.58 % instead of 51.67 % in the audit example). U is now counted as
  T. Since the menu gate and `core.api` already hand RNA to this module
  in the DNA alphabet, only direct module calls were affected.
- Tests: 3 new (`ModuleInputTests`). Suite 194 → 197. Golden master
  unchanged.

With this package every D14 sub-finding (F1–F12) and D13 is corrected.
A re-run of the read-only probe from the investigation shows each module
either accepting a type or refusing it with a message; no module returns
a silent result for a foreign type any more.

## Phase 6 – Tier 2 functional corrections (2026-10-31)

Tests are in `tests/test_phase6_tier2.py`.

### D6, S7 Programmatic local BLAST with a database path (`analyses/blast.py`)

- Problem: `run_blast_analysis(engine="local", local_database_path=...)`
  tested the path with `Path.exists()`. A BLAST database is a path
  *prefix* (`.../synapse_test` for `synapse_test.nin`, ...), so the check
  always failed for a valid database. On that branch the database type
  was also guessed from the catalog name, and title, date, version and
  size stayed `None` (S7).
- Now: the prefix is scanned and verified by `LocalBlastDatabaseManager`
  (rooted at the prefix's folder), exactly like a managed database:
  type and volumes from the files, title/date/sequence count from
  `blastdbcmd -info`. A missing prefix raises a clear error, and a
  program/database type mismatch is refused as on the managed path.
- Validation: external tool execution — the bundled `data/local_test`
  database addressed by its prefix gives the two hits of the audit's
  first real run (coverage 100 % and 98.04 %) with full provenance.
- Tests: 3 new (`LocalBlastDatabasePathTests`). Suite 197 → 200.

### D7 The stale BLAST artifact is labeled as legacy (`blast_results/README.md`)

- Problem: `blast_results/TEST_8_REFERENCE_blast.json` still carries the
  pre-P0-1 values (100 % coverage for a 32 bp HSP on a 99 bp query, E =
  7.45 without a band) in an obsolete key schema, and
  `docs/AUDIT_STATUS.md` cited it as the evidence for the P0-1 fix.
- Now: the file is kept unchanged as the record of a real NCBI run (RID
  `CD3MMW95016`) and of the defect; `blast_results/README.md` states that
  it is not a valid result, what is wrong in it, and where the corrected
  values are tested. `docs/AUDIT_STATUS.md` P0-1 cites the test, not the
  file. Deleting the file was not done: it removes the only record of
  that remote run and needs the owner's decision.

### D9 The guide's own site is not an off-target (`analyses/crispr.py`)

- Problem: the local specificity screen counted the guide's intended site
  as an exact off-target. Supplying the designed sequence as the local
  reference - what the menu's "Off-Target Analysis" invites - gave every
  guide `exact_matches = 1` and a specificity component of 60 instead of
  100, which entered the ranking.
- Now: an exact hit on the same strand at the candidate's own guide
  coordinates is marked `on_target` (classification "on-target"). It
  stays in the hit list, is reported as `on_target_matches`, and is not
  counted in `exact_matches`, the mismatch classes, `total_hits`,
  `off_target_hits`, the specificity score or the "local matches" flag.
  At most one hit per candidate is the on-target site, so an identical
  copy in a second reference still counts as an off-target.
- Limitation: the rule compares coordinates as given; it identifies the
  own site when the reference uses the coordinates of the designed
  sequence (the designed sequence itself, or its parent).
- Tests: 3 new (`CrisprOnTargetTests`, both strands). Suite 200 → 203.

### D8 CRISPR guides report the modeled cut site (`analyses/crispr.py`)

- Problem: no guide carried a cut position; the only cleavage information
  was a free-text `cut_model` string.
- Now: each system defines where it cuts, counted along the guide in
  targeting orientation, with its source: SpCas9 and SpCas9-NG blunt
  after guide position 17 (3 nt 5' of the PAM; Jinek et al. 2012), SaCas9
  blunt after position 18 (3 nt 5' of the PAM, modeled as for SpCas9;
  Ran et al. 2015), Cas12a staggered after position 18 on the PAM strand
  and 23 on the target strand, a 5-nt 5' overhang (Zetsche et al. 2015).
  Every `GuideCandidate` has `cut_site`, `cut_site_complementary` (the
  break lies between forward positions p and p+1, 1-based) and
  `cut_type`; they appear in the candidate report, the system details
  (with the citation) and the JSON export.
- Scope: this is the canonical cleavage geometry, not a prediction of
  cleavage efficiency or of the repair outcome.
- Validation: literature (PMIDs 22745249, 25830891, 30166441, 26422227,
  checked in PubMed); tests check on both strands that the guide prefix
  up to the cut lies 5' of the break in the sequence.
- Tests: 3 new (`CrisprCutSiteTests`). Suite 203 → 206.

### D10, D11 CRISPR configuration fields act; exports stay in the project (`analyses/crispr.py`)

- D10 problem: `include_low_quality`, `rank_by_gc`, `editing_window_start`,
  `editing_window_end` and `base_editor` in `DesignConfig` were never
  read.
- Now: `include_low_quality=False` removes candidates with FAIL flags and
  says how many in `warnings`; `rank_by_gc=False` takes the GC component
  out of the rank weights (the remaining weights are re-normalized; GC is
  still reported); the configured editing window is used by the menu's
  base-editing analysis and by `find_base_editing_candidates`, which also
  takes `base_editor` when no editor is passed (`config_editing_window`
  validates `1 <= start <= end`). The defaults equal the previous
  behavior (window 4–8 = the CBE/ABE definitions; GC ranked; low-quality
  candidates included), so default results are unchanged.
- D11 problem: relative CRISPR export paths resolved to the user's
  Desktop (here a OneDrive-synced folder).
- Now: they resolve to `exports/crispr/` in the project (ignored by Git,
  next to the report exports), or to `$SYNAPSE_EXPORT_DIR/crispr`.
- Also fixed: `analyses/crispr.py` did not import `os`, needed by the new
  export-directory code (found by the new test before commit).
- Tests: 4 new (`CrisprConfigurationTests`). Suite 206 → 210.

### D19 The pipeline's ORF runner takes the menu's search settings (`core/pipeline.py`)

- Problem: `SynapsePipeline.run_orf` exposed none of `start_codon_policy`,
  `nested`, `ranking` or `genetic_code`, so a programmatic pipeline run
  could not reproduce an interactive ORF search. `PipelineContext.completed`
  started with six keys and grew with whatever was stored, so `status()`
  had a history-dependent shape while `print_status` listed twelve.
- Now: `run_orf` takes the four settings (None = module default) and
  passes them to the ORF module; its results equal `core.api.find_orfs`
  for the same settings. `completed` always has the twelve modules of
  `STATUS_LABELS` (one table, also used by `print_status`).
- Tests: 2 new (`PipelineOrfRunnerTests`). Suite 210 → 212.

### D15 The unimplemented planner is labeled as such (`main.py`, `ui/menu.py`)

- Problem: the main menu advertised `[11] AI Experiment Planner` and only
  printed "coming soon"; `docs/AUDIT_STATUS.md` lists it as P3-1
  DEFERRED.
- Now: the entry reads `[11] AI Experiment Planner (not implemented)` and
  selecting it states that nothing is run or stored and that the planner
  is deferred (P3-1). The key is kept so the menu numbering is unchanged.
  Implementing the planner is not part of this correction.
- Tests: 1 new (`PlannerPlaceholderTests`, drives `main.py`). Suite 212 → 213.

### S15, S16 Two variant-module defects (`analyses/variant_analysis.py`)

- S15 (confirmed by execution): `annotate_from_mutation_analysis` always
  failed - it imported a top-level `mutation` module that does not exist
  (`analyses.mutation` is the module) and would then have passed
  `reference_cds`/`variant_cds`, which `compare_sequences` does not
  accept. It now imports `analyses.mutation` and calls `compare_sequences`
  with its real parameters. pyright's two diagnostics on that call are
  gone (baseline 134 → 132).
- S16: `left_normalize_variant` added `REF_MISMATCH` through a shallow copy
  of the annotation, so the flag also appeared in the *input* variant's
  existing flag list. It now uses the copy-on-write helper introduced for
  D12.
- Tests: 2 new (`VariantModuleFindingTests`). Suite 213 → 215.

### S17 The cloning digest cuts where the enzyme cuts (`analyses/cloning.py`)

- Problem: `_selected_digest` (used for the vector and insert digests in
  the cloning strategy) used each site's recognition start as the cut
  position, so fragment lengths were off by the enzyme's top-strand cut
  offset minus one (EcoRI agreed only because its offset is 1). Audit
  example: cuts 51/107 and fragments 51/56/55 bp instead of 51/111 and
  51/60/51 bp.
- Now: the boundary is recognition start − 1 + `cut_top` (modulo the
  length for circular molecules), the convention of
  `restriction.calculate_cut_positions`; the recognition start is kept
  in each cut record.
- Validation: for all 231 enzyme pairs of the database, linear and
  circular, the cloning digest's fragment lengths equal those of
  `restriction.perform_digest`, whose circular digest agrees with
  Biopython `Bio.Restriction` (D3).
- Tests: 2 new (`CloningDigestCutTests`). Suite 215 → 217.

### S19 Start-codon classification and heuristic scores follow the genetic code (`analyses/orf.py`)

- Problem: `classify_start_codon`, `calculate_cds_quality` and
  `calculate_orf_confidence` used the fixed ATG/GTG/TTG set. A CTG ORF
  found with `table_starts` was classified "None", got no "Alternative
  start codon" reason and a CDS-quality score of 60 instead of 70 for an
  otherwise identical GTG ORF; table-2 starts (ATA, ATT) were treated the
  same way.
- Now: all three use `initiation_codons(genetic_code)` (the D1 rule);
  `classify_start_codon` takes an optional `genetic_code` (default 1).
  Results for ATG, GTG and TTG ORFs are unchanged; the golden master is
  unchanged.
- Tests: 2 new (`StartCodonClassificationTests`). Suite 217 → 219.

### S23, S3, S4, S5 CRISPR seed orientation, RTT coverage and sourced thresholds (`analyses/crispr.py`)

Performed by Opus 5 on a restricted assignment, after a model safeguard
interrupted the Opus 5.5 session at the start of this investigation.

- **S23 (confirmed by execution, corrected).** The stored `seed_sequence`
  took the guide's 3' end for every system. That is the PAM-proximal seed
  only for 3'-PAM systems; for Cas12a (5' PAM) it reported the PAM-distal
  end. Reproduced: guide `CGTTACGCCCGTGGACAGAATTA` stored `CAGAATTA`
  instead of `CGTTACGC`. The local specificity screen already branched on
  `pam_side`, so the stored seed and the screened seed disagreed. Both now
  use the new `pam_proximal_seed()`, which is the single source of the
  rule. Tests check for all four systems that the seed sits directly
  against the PAM.
- **S3 (confirmed, corrected).** `edit_in_rtt` only compared the RTT length
  with the requested length and therefore could not state whether the RTT
  covers the edit site. It is now a coordinate check against the clipped
  RTT window. The emission filter already tested the RTT length separately,
  so this tightens the filter without loosening any existing condition. The
  report line is now "RTT spans edit site". Remaining limitation, stated in
  every candidate's notes and in METHODS: the RTT is the reference template
  and does not encode the alternative allele, so it is not a ready-to-use
  pegRNA RTT. Renaming the exported field belongs to D18.
- **S4 (resolved by documentation).** No threshold was changed. The GC
  range, homopolymer limit and local-mismatch limit are now labeled Synapse
  screening conventions; the seed length cites the PAM-proximal seed concept
  (Sternberg et al. 2014) while labeling 8 nt as a Synapse choice; the
  base-editing windows cite Komor et al. 2016 (CBE) and Gaudelli et al. 2017
  (ABE), with the `+` variants labeled as Synapse windows.
- **S5 (partly confirmed, addressed).** The audit's concern that the poly-T
  flag is position-independent is not a defect: Pol III terminates wherever
  it meets the T-tract. The real gap was that the motif was neither named
  nor sourced. The flag now reads "WARN: poly-T run (possible U6/Pol III
  terminator)" and the threshold is a named constant. METHODS records that
  four T is a conservative screening threshold, not a prediction, citing
  Verosloff et al. 2021: four uracils terminate at Pol III type 2
  promoters, whereas type 3 promoters such as U6 required longer tracts or
  upstream RNA secondary structure.
- Validation: literature verified in PubMed (PMIDs 24476820, 27096365,
  29160308, 31634902, 33811918). No external tool execution was required.
- Tests: 9 new (`CrisprSeedOrientationTests`, `PrimeEditingRttTests`,
  `CrisprThresholdDocumentationTests`). Suite 219 -> 228; pyright unchanged;
  golden master unchanged.

## Phase 6 - Architecture (2026-10-31)

### D16 (part 1) One definition per name in the primer module (`analyses/primer.py`)

- Problem: `calculate_tm_nearest_neighbor` was defined twice at module
  level. The first definition was the V3 model; it was then rebound to
  `_calculate_tm_nearest_neighbor_v3_raw` and the public name was replaced
  by a backward-compatible adapter. The behavior was correct, but the
  definition a reader finds first was not the one that runs, and pyright
  reported the shadowing as `reportRedeclaration`.
- Now: the model is defined directly under its private name and the alias
  assignment is gone, so the public name has exactly one definition - the
  adapter - and the module has no duplicated module-level name at all.
- Equivalence: all 129 module-level callables of `analyses.primer` were
  fingerprinted by the SHA-256 of their compiled bytecode before and after.
  No name was added or removed and no entry changed; both Tm names keep
  exactly the code object they had. Module-level lookups resolve at call
  time, so every existing call site reaches the adapter as before.
- pyright: the `reportRedeclaration` diagnostic is gone (134 -> 131 against
  the `ab4bd2a` baseline, together with the two S15 diagnostics).
- Test changed deliberately:
  `test_primer_has_single_live_definition_per_name` asserted
  `duplicates == {"calculate_tm_nearest_neighbor"}`, i.e. it recorded the
  duplication as intended. It now asserts `set()`, which is stricter, and a
  new test checks that the adapter's legacy parameter aliases reach the same
  model as the v3 names.
- **Not done, deliberately:** the live `_v5`-`_v11` ranking wrapper chain is
  unchanged. Flattening it alters executed code and requires the independent
  review recorded as open validation item 7 in `docs/AUDIT_STATUS.md`. D16 is
  therefore partially corrected.

### S8, S10, S21, S22 Four silent-failure paths now report what happened

- **S21** (`fileio/genbank.py`): the ORIGIN reader keeps only letters, which
  is correct for the line numbers and block spacing of the format but also
  discarded any other symbol without a word. Symbols that are neither
  letters, digits nor whitespace are now collected and reported as a record
  note, so a sequence shorter than the file suggests is visible.
- **S10** (`fileio/genbank.py`): locations in forms the reader does not
  parse - `order()`, remote accession references, between-base sites - keep
  only their raw string and yield no interval. The record now carries a note
  naming how many such locations there are and which feature keys they
  belong to. A record whose locations all parse gets no note; the bundled
  fixture is unchanged.
- **S8** (`core/report.py`): `to_jsonable` caught every exception from a
  payload's `to_dict()` and fell through to `repr()`, so an unparseable
  string could reach a report with nothing to explain it. The failure is now
  logged as a warning with the traceback, and the generic repr fallback is
  logged at debug level. The serialized output is unchanged.
- **S22** (`analyses/variant_analysis.py`): the supplied reference was never
  checked, so a sequence in another alphabet was compared against REF
  alleles base by base and produced mismatches that looked like data
  problems. `validate_reference_sequence()` and `validate_contig_references()`
  now require the IUPAC DNA alphabet and are applied by `normalize_variants`
  and `check_reference_alleles`. Ambiguity codes stay valid, since real
  assemblies contain N and occasionally other IUPAC codes; U is refused,
  which matches the VCF specification's ACGTN alleles.
- Tests: 7 new (`GenBankReaderNoteTests`, `ReportSerializationTests`,
  `VariantReferenceAlphabetTests`). Suite 229 -> 236; pyright unchanged;
  golden master unchanged.

### S1, S2, S9, S11, S12 Five suspected findings assessed

- **S1 is not a defect.** The Kozak context of a minus-strand ORF is
  computed in the reverse-complement frame, which is the biologically
  correct frame: that frame *is* the minus strand read 5'->3'. Reproduced:
  a transcript carrying `GCCACC` upstream of its start, placed on the minus
  strand of the parent, is reported with context `GCCACCATGG` and
  classification "Strong". The result carries no coordinate, so nothing can
  be mapped into the wrong frame. The frame is now documented and locked by
  a test; no code was changed.
- **S2 is not a defect.** `_keep_longest_per_stop` collapses minus-strand
  ORFs correctly: two in-frame starts sharing one stop get the same frame
  and the same 3'-end key, and only the longest survives with
  `nested=False`. The audit's concern about ORFs "found from different
  rotations" does not arise, because one call uses one rotation of the
  reverse complement. Locked by a test; no code was changed.
- **S9 is confirmed and the note now quantifies it.** Peak memory of the
  global alignment was measured with `tracemalloc` at 300, 600 and 1,200 nt
  and scales quadratically at roughly 45 bytes per matrix cell. The note
  therefore states about 180 MB at 2,000 nt and about 4.5 GB at 10,000 nt,
  instead of only saying that memory "grows". The algorithm is unchanged;
  a banded or linear-space alignment remains open as P2-6.
- **S12 is resolved.** pyright was not run during the audit. It is now run
  against the 134-diagnostic baseline of `ab4bd2a` for every work package.
  Note for future sessions: `pyright.exe` is blocked by a Windows
  application control policy (WinError 4551), so it must be invoked as
  `python -m pyright --pythonversion 3.10`.
- **S11 remains unverified, by choice.** The literature module's behavior
  against live providers (rate limits, partial outages, cross-source
  duplicate DOIs) was not tested, because that means sending queries to
  external services; the resource policy requires approval for that. The
  offline fixtures continue to cover parsing, deduplication and retry
  handling. This is recorded as a limitation, not as a defect.

### D16 (part 2) The superseded primer generations are gone (`analyses/primer.py`)

- Investigation first, as the audit required. Three independent checks
  established what the `_vN_` chain actually does: an AST pass over every
  module-level binding in source order; a project-wide reference scan; and a
  `sys.settrace` run of a real `design_primers_professional()` call. The AST
  pass corrected an assumption carried in the audit - `design_primers_professional`
  is **not** `_v7_design_primers_professional` but a plain definition near the
  end of the module that overrides the v7 alias.
- Result: each public name was rebound up to five times in sequence
  (`rank_professional_primer_pairs` ran through v5, v8, v9, v10, v11), so a
  reader had to scan the whole module to learn which implementation runs.
  Only the last binding was ever reachable.
- Removed: the 7 superseded alias statements and, by repeating the reference
  scan to a fixpoint, the 24 private definitions that only those aliases kept
  reachable. `analyses/primer.py` goes from 6,434 to 5,108 lines.
- Safety: there are no import-time calls in the module (checked by AST: the
  only non-declarative top-level nodes are a dict update and two
  `if __name__ == "__main__":` guards) and no dynamic lookup of module
  attributes, so an intermediate binding could never have been observed.
- Equivalence: the SHA-256 of the compiled bytecode of every module-level
  callable was compared before and after. 129 -> 105 callables, **nothing
  added, nothing changed**, 24 removed - every surviving function keeps its
  exact code object. The four live bindings are unchanged.
- pyright: 134 -> 112 against the `ab4bd2a` baseline. The 22 removed
  diagnostics all belonged to the deleted code; no new diagnostic appeared.
- A reading guide above the ranking section now states which public name is
  bound to which implementation and that an older `_vN_` prefix on a
  surviving helper is history, not a version switch.
- Not done: merging the surviving generations into single functions. That
  would change executed code and has no reachability argument behind it, so
  it is out of scope.

### P1-9 The RNA codon view uses the selected ORF's genetic code (`analyses/rna.py`)

- Problem, the last remaining gap of P1-9: the three-frame codon preview
  always read NCBI table 1. With a table-2 ORF selected, the preview stopped
  at the first UGA and showed `M`, while the ORF's own translation correctly
  read UGA as Trp and gave `MWWWWWK`. The same screen therefore showed two
  different answers for one sequence.
- Now: `rna_codon_table()` returns the table for a code,
  `translate_rna_frame()` and `analyze_rna_three_frames()` take
  `genetic_code`, and `analyze_rna()` records it in the result. The chosen
  table reaches them from the selected ORF: `orf_genetic_code()` reads it
  from the ORF, `analyze_rna_from_pipeline()` and `SynapsePipeline.run_rna()`
  apply it unless the caller passes another one, and `core.api.analyze_rna`
  exposes it as a parameter.
- Verified agreement: for a table-2 ORF the +1 frame preview now equals
  `orf.translate_orf_sequence(orf)` exactly, which is what the finding asked
  for. Without a selected ORF the default stays table 1, so existing results
  and the golden master are unchanged.
- Tests: 5 new (`RnaGeneticCodeTests`). Suite 240 -> 245.

### S18 The cloning planner stops repeating identical site scans (`analyses/restriction.py`)

- Measured cause, by profiling a 1.3 kb vector and a 0.3 kb insert with 8
  enzymes: 1,240 evaluated strategies produced **5,002** calls to
  `find_recognition_site_matches` for only **16** distinct combinations of
  sequence, recognition sequence and topology. The scan itself was 98 % of
  the runtime. The planner's results were never wrong; it redid the same
  work about 312 times over.
- Now: the scan is a private `_scan_recognition_sites()` behind an
  `lru_cache`, keyed by exactly the three values it depends on.
  `find_recognition_site_matches` rebuilds fresh dicts from the cached
  tuple, so a caller that edits a returned match cannot reach into the
  cache.
- Measured effect (same input): 8 enzymes 11.8 s -> 0.24 s; 10 enzymes,
  which the audit timed at 100-110 s, 0.57 s; all 22 enzymes, which the
  audit timed at over two minutes, 18.7 s. The remaining time is the
  strategy enumeration itself, not repeated scanning.
- No scientific behavior changed: the scan is a pure function of its three
  inputs, the full suite including the S17 comparison across all 231 enzyme
  pairs passes unchanged, and the golden master is untouched.
- Tests: 3 new (`CloningPlannerScanReuseTests`), covering the scan count,
  isolation from caller mutation, and identical cold/warm results for both
  topologies. Suite 245 -> 248.

### S14 The API covers the remaining analyses (`core/api.py`)

- Problem: `core/api.py` offered ten entry points but none for cloning,
  variant filtering, mutation reporting, BLAST or literature, so those
  workflows were reachable only through the terminal menus.
- Added: `plan_cloning`, `filter_variants`, `mutation_report`,
  `blast_search` and `literature_search`. Each returns an `AnalysisResult`
  with its parameters, warnings and limitations, like the existing ten.
- Privacy is not weakened by the move off the menu. The BLAST menu asks
  before a remote search sends the query to NCBI; the API cannot ask, so
  `blast_search(engine="remote")` raises unless the caller passes
  `send_to_ncbi=True`, and `literature_search` likewise requires
  `send_to_providers=True`. A script therefore cannot transmit sequences
  or queries unintentionally. The local BLAST engine needs no confirmation
  because it sends nothing.
- `plan_cloning` accepts DNA only, since restriction enzymes cut DNA, and
  takes each record's topology unless told otherwise.
- Tests: 8 new (`ApiCompletenessTests`), including both privacy gates and a
  local BLAST+ run that sends nothing. Suite 248 -> 256.
- Still open for P2-2: the analysis modules continue to own their terminal
  menus. The API now covers every analysis, which is the part S14 named;
  moving the menus out of the modules belongs with the Stage F interface
  work.

### S6, P1-7 An NCBI search can be resumed from its request ID (`analyses/blast.py`, `core/api.py`)

- Problem: every run already recorded its RID and the client had the
  primitives (`submit`, `status`, `fetch_xml`), but nothing used them to
  pick a search up again, and no menu entry offered it. A long search that
  outlived its session had to be resubmitted.
- Now: `NCBIBlastClient.await_ready()` holds the polling loop that `run()`
  always had, and `NCBIBlastClient.resume(rid)` waits for an existing
  request and parses its result. `run()` now submits and then calls
  `resume()`, so both paths share one implementation. Menu entry
  `[14] Resume Search by Request ID (RID)` and `core.api.blast_resume()`
  expose it.
- Nothing is invented for a resumed run. `describe_run_from_xml()` rebuilds
  the query and parameters from what the XML states - definition, length,
  program, database - and leaves the query sequence empty, because BLAST
  output does not contain it. `BlastQuery.reported_length` carries the
  stated length so coverage still works; when the XML states no length,
  the query says so instead of guessing. A resumed run is marked as such
  in its warnings, since the original submission time is unknown.
- Privacy: resuming sends only the request ID, never the sequence again.
  The menu asks before contacting NCBI and `blast_resume` requires
  `contact_ncbi=True`, consistent with the other external entry points.
- An expired request ID reports NCBI's `UNKNOWN` status with an
  explanation, because NCBI retains an RID only for a limited period.
- Tests: 9 new (`BlastResumeTests`), driven by a stubbed client that
  asserts `submit()` is never called, so **no request leaves the machine**.
  Suite 256 -> 265.
- Not verified: a real resume against NCBI. That needs a live submission,
  which requires the owner's approval, so it remains an untested path
  against the real service.

### S13 The session runs from `main()`, not at import (`main.py`)

- Problem: `main.py` held its control flow at module level - a call to
  `choose_sequence_file()` followed by a `while True:` loop - so importing
  the module started the application. Tests could only drive it in a
  subprocess, and any accidental `import main` would have blocked on input.
- Now: that flow is the body of `main()`, called from an
  `if __name__ == "__main__":` guard. The 650 lines were indented
  unchanged; `git diff -w` reduces the whole change to the 12 added lines
  of the function header, its docstring and the guard, which is the
  equivalence evidence for a move of this size.
- The helper functions and module constants stay at module level, so
  `report_error`, `choose_sequence_file`, `run_protein_record` and the rest
  remain importable and testable on their own.
- Tests: 2 new (`MainEntryPointTests`): a fresh interpreter imports `main`
  with no stdin and the session does not start, and the module has no
  top-level loop and exactly one `__main__` guard calling `main()`. The
  existing smoke tests still drive `python main.py` end to end.
  Suite 265 -> 267.

## Phase 6 - Stage F: terminal UI standardization (2026-11-01)

### D17 One visual language for every terminal interface

A mechanical inventory of `analyses/`, `ui/`, `core/` and `main.py`
measured the problem before anything was changed: **11 separator widths**
(42, 45, 58, 60, 65, 68, 70, 72, 78, 88, 90) across 389 literal sites,
**7 meanings for `[0]`**, and three spellings of the error prefix. BLAST
framed its screens at 90 columns, literature at 88, CRISPR at 78.

- **`ui/style.py`** is the new single source of the visual language:
  rules, headers, sections, aligned fields, menu entries, navigation
  entries, input prompts and message prefixes. It holds no scientific
  logic and imports nothing from `analyses/`, `core/`, `data/` or
  `fileio/`.
- **One width.** `RULE_WIDTH` is 70 - the width already most common, and
  it fits an 80-column terminal. 382 literals in 16 files now call
  `major_rule()` or `minor_rule()`; the inventory reports zero remaining
  distinct widths. `blast._line()` and `crispr._line()` stay as local
  shorthands but default to the shared width instead of 90 and 78.
- **One meaning for `[0]`: leave the screen.** `Back` in a menu, `Cancel`
  in a dialog that would otherwise commit, `Exit` in the main menu only.
  Four entries were corrected: `orf.py` "Apply and return" (settings are
  applied as they are changed, so it was always a Back), `literature.py`
  "Back / Exit", `main.py` "Choose another sequence" (now "Cancel" with
  its effect named), and `mutation.py` "Enter sequence manually", which
  was a choice rather than an exit and moved to `[M]`.
- **One prefix per message category:** `ERROR:`, `WARNING:`, `NOTE:`,
  `OK:`. Eight `Error:` and one `Note:` variant were corrected. Scientific
  wording stays module-specific; result flags such as `WARN: poly-T run`
  are data, not interface messages, and are untouched.
- **Two deliberate exceptions, documented in place and in the standard:**
  the ORF candidate table draws a rule that spans its own table and uses
  the standard width as a floor; `core/pipeline.py` prints terminal output
  from inside `core/` and takes `ui.style` through function-local imports
  so `core` keeps no module-level dependency on `ui`.
- **One known deviation, recorded rather than hidden:**
  `mutation._prompt_sequence` now offers no `[0]`, because a real cancel
  needs a cancel contract and four call sites changed - a control-flow
  change, not a presentation one.
- `docs/TERMINAL_UI_STANDARD.md` documents the standard with the measured
  before/after. `tests/test_terminal_ui.py` enforces it across the whole
  project: no separator literal outside `ui/style.py`, `[0]` restricted to
  the three navigation verbs, `Exit` only in the main menu, no prefix
  variants, and the layering of `ui.style` and `core`. The checks were
  verified against deliberate violations before being relied on.
- No wording of a scientific statement, no menu structure and no
  analytical output changed. Suite 267 -> 277; golden master unchanged;
  pyright 112 with no new diagnostics.

### Stage F: the BLAST menu reads in order again

The RID resume entry added for P1-7 had been given the next free number,
`[14]`, inside the SEARCH block, so that block read 1, 2, 3, 14. Stage F
is where a coherent renumbering belongs, so the whole menu was reordered
to follow the work rather than the order features were added: run,
resume, configure, manage; then result analysis; then output; then
system. Numbers now run 1-14 with `[0] Back`.

Both the printed labels and the dispatch branches were rewritten from one
mapping, and each branch was checked against its new label. The rewrite
was scoped to `blast_menu` by its AST bounds, because `analyses/blast.py`
contains several other numbered sub-dialogs that must not shift.

`tests/test_terminal_ui.py` now also checks that a menu's numbers are
contiguous and that every entry has a branch and every branch an entry,
for both the BLAST menu and the main menu.

### S24 The XML parsing exposure is measured, not assumed

The finding was recorded as a suspicion. It is now measured on the
interpreter in use (Python 3.14.6, Expat 2.8.1), and the result is mixed
rather than a clean pass or fail:

- **XXE is not reachable.** `xml.etree.ElementTree` does not resolve
  external entities: an entity pointing at a local file raises
  `ParseError: undefined entity`, and an external DTD is ignored rather
  than fetched.
- **Entity expansion is reachable but capped.** A billion-laughs document
  does expand - 1 MB in 0.02 s - until Expat's input-amplification guard
  stops it (`limit on input amplification factor ... breached`, 10 MB
  blocked in 0.09 s). The practical cost is about a megabyte of work, not
  a denial of service.
- **The residual risk is version-dependent.** That guard arrived in Expat
  2.4.0. Synapse targets Python 3.10, and a 3.10 build bundling an older
  Expat would expand without a cap. The test skips the amplification check
  on such a build rather than reporting a pass it did not observe.

`defusedxml` is still not justified as a runtime dependency under the
standard-library-first policy, and the generic advice to adopt it does not
match what was measured here. If the project ever has to run on a build
without the guard, setting the limits explicitly through `xml.parsers.expat`
is the dependency-free hardening.

Tests: 4 new (`XmlParsingExposureTests`).

## Phase 6 - Stage G: scientific output standardization (2026-11-02)

### D18 Heuristic scores are no longer called "confidence" (`analyses/orf.py`, `analyses/protein.py`, `analyses/cloning.py`, `analyses/primer.py`)

The ORF display was relabelled "Heuristic score" in Stage 1, but the data
model and every export still called the value `confidence`, a word that
implies a probability. The value is a sum of fixed points (complete ORF
+25, ATG +20, length bands, translation safety, ambiguous codons, Kozak,
Shine-Dalgarno), clipped to 0-100. It ranks candidates; it does not
estimate how likely a candidate is to be coding.

| Where | Before | After |
|---|---|---|
| ORF record | `orf["confidence"] = {"score", "classification"}` | `orf["heuristic_score"] = {"score", "classification"}` |
| Scoring function | `calculate_orf_confidence` | `calculate_orf_heuristic_score`; the former name stays as an alias |
| `compare_orf_candidates` rows | `confidence_score`, `confidence_classification` | `heuristic_score`, `heuristic_classification` |
| `calculate_orf_feature_summary` | `high_confidence`, `moderate_confidence`, `low_confidence` | `high_heuristic_score`, `moderate_heuristic_score`, `low_heuristic_score` |
| ORF interpretation note | "Confidence and CDS quality scores are heuristic ranking systems" | "Heuristic and CDS quality scores are ranking systems" |
| Signal-peptide screen | `"confidence": "low"/"moderate"`, printed as "Confidence:" / "Confidence label:" | `"evidence": "low"/"moderate"`, printed as "Heuristic evidence:" |
| Cloning plan score component | `metadata_confidence` | `compatibility_metadata` (it awards 10 points when the end compatibility is known) |
| Primer warning when every pair is MARGINAL or POOR | "No high-confidence primer pair was found ..." | "Every reported primer pair is rated MARGINAL or POOR ..." |

Values, bands, thresholds and ranking are unchanged. Only names changed.

**Golden master (§11.7).** `tests/golden/orf.json` and
`tests/golden/orf_legacy_baseline.json` stored the score under the key
`confidence` (6 and 23 occurrences). The key is now `heuristic_score`.
The snapshots were rewritten by renaming that key and nothing else; with
the key mapped back, both files are identical to the previous version,
and the unmodified test then confirmed that the code produces exactly the
renamed files. The large line count of the diff comes from the snapshot
writer's sorted keys: the renamed key moves from between `cds_quality` and
`frame` to between `frame` and `length`.

| Snapshot | SHA-256 before | Keys renamed |
|---|---|---|
| `orf.json` | `e5fe5dd9...5037a672` | 6 |
| `orf_legacy_baseline.json` | `67fef456...56d77627` | 23 |

**Compatibility.** Code that reads `orf["confidence"]` or the old
comparison and summary keys must use the new names; no alias is kept in
the data, because an alias would put the misleading word back into every
export. The function name is aliased because an alias there does not
reach any output.

Unchanged on purpose: docstrings that say a value is *not* a statistical
confidence (BLAST evidence profile, literature relevance score, protein
sequence-applicability label) use the word correctly.

Tests: 5 new (`HeuristicScoreVocabularyTests`), each confirmed to fail on
the previous code. One of them checks every string literal used as a
dictionary key or subscript in `analyses/` for the word.

### D18 (part 2) The sequence-type classifier reports an evidence level (`core/alphabets.py`, `core/models.py`, `core/report.py`, `fileio/loader.py`, `main.py`, `analyses/blast.py`)

The Phase 6 classifier (D13/D14) called its CONFIRMED / PROBABLE /
AMBIGUOUS / INVALID outcome a "confidence". Those four values are
categories of evidence about which letters occur; none of them is a
probability. The word is now reserved for statistical confidence, which
Synapse does not compute.

| Where | Before | After |
|---|---|---|
| `SequenceClassification` field | `confidence` | `evidence` |
| `SequenceClassification.as_dict()` | `"confidence"` | `"evidence"` |
| Level constant | `CONFIDENCE_LEVELS` | `EVIDENCE_LEVELS` |
| `SequenceRecord` attribute | `type_confidence` | `type_evidence` |
| Provenance and report input section | `sequence_type_confidence` | `sequence_type_evidence` |

The values and the classification rules are unchanged; no golden snapshot
contained these names. The earlier entries "Central sequence alphabets and type
classifier" and "D14-F3, F11 Declared molecule type is read, kept and
reported" describe the names as they were introduced and are left as the
record of that change.

Tests: the existing sequence-type tests now read the new names, and the
source check from part 1 covers `core/`, `fileio/`, `ui/` and `main.py`
as well as `analyses/`.

### Stage G One notation for units and for primer coordinates (`analyses/primer.py`, `analyses/blast.py`, `analyses/rna.py`)

An inventory of every f-string literal in `analyses/`, `core/`, `fileio/`,
`ui/` and `main.py` found one dominant notation per unit and a few
deviations:

| Quantity | Dominant notation | Deviations | Now |
|---|---|---|---|
| Percentage | `45.20%` (59 sites) | BLAST identity and coverage `45.20 %` (2) | `45.20%` |
| Temperature | `60.00 °C` (10 sites) | primer `ΔG @25°C` (3) | `ΔG at 25 °C` |
| Entropy | - | `cal/mol/K` (2) and `cal/(mol*K)` (3) | `cal/(mol·K)` |
| Molecular mass | `Da` (protein, 5; DNA as `kDa`) | RNA `g/mol` (2) | `Da` |

Values and precision are unchanged; only the unit text changed. Da and
g/mol are numerically identical here, so the RNA figures are the same
numbers under the unit the other modules use.

**Primer coordinates.** `docs/METHODS.md` states that user-facing
positions are 1-based and inclusive, and the ORF, CRISPR, variant and
main menus follow it. The primer menu did not: it printed the internal
0-based half-open `start-end` unlabelled, so a primer covering bases 117
to 135 was shown as `116-135`, and it asked for the target region as
0-based half-open. It now prints `117-135 (1-based)` for primers and the
target region and asks for 1-based inclusive positions, which it converts
to the half-open interval the design functions take. Out-of-range input
is refused with the 1-based rule (`1 <= start <= end <= length`). The
design functions, `core.api` and exported results keep their 0-based
half-open fields, so no programmatic caller changes.

Not changed, recorded instead: in `analyses/restriction.py` the key
`site_position` is 1-based in the restriction map and 0-based in digest
cut records (which also carry explicit `site_position_1based`). The key
is part of `tests/golden/restriction.json` and of the results compared
with Biopython on pUC19; making it uniform changes stored values, which
is an interface change of its own and not a notation fix.

Tests: 3 new (`OutputNotationTests`): a source check of every f-string
literal for the four unit rules, the 1-based primer display, and the
1-based target input with six cases including out-of-range values. Each
fails on the previous code.

### Stage G (6.4) Warnings a result carries reach the exported summary (`core/models.py`)

Each stored analysis is exported with its parameters, warnings and
limitations. The JSON report holds the full payload as well; the Markdown
report prints only those three lists. A measurement over every offline
`core.api` analysis showed that two payloads carry warnings of their own
that never reached the list: a primer design on a 600-nt test sequence
had three (for example "Specificity was not genome-wide screened; no
reference database was supplied."), and CRISPR designs carry theirs the
same way. Through `core.api` the Markdown report therefore showed none of
them. The interactive primer path copied them by hand, so the two paths
disagreed.

`AnalysisResult` now adds a payload's top-level `warnings` to its own
list when it is created, without duplicates. Every path that stores a
result uses it, so the API and the menus agree, and the payload and its
`data_sha256` are unchanged.

Measured at the same time: results stored from the menus for sequence
overview, restriction and cloning recorded no parameters, while the API
path records topology, enzymes and the other inputs, and neither path
recorded which version of the enzyme definitions was used. Restriction
results from the menu now record `enzymes` and `topology`, named
explicitly from the analysis result, and restriction and cloning results
on both paths record `enzyme_database_version`
(`data/restriction_enzymes.py`, currently `synapse-2026-10-29`). The
menu-path sequence overview and cloning analysis still record no input
parameters; that remains follow-up work, because naming them means
choosing per module which inputs matter rather than copying the payload.

Tests: 3 new (`ExportedWarningsTests`), each failing on the previous
code.

## Phase 6 - Stage H: language and documentation (2026-11-02)

### Editorial review, part 1: spelling, recurring messages, tone

The review began with a mechanical inventory of every string literal in
`analyses/`, `core/`, `fileio/`, `ui/` and `main.py` (docstrings counted
separately) and of every Markdown document, against five questions:
German text in the interface (none found), mixed spelling, casual tone,
claims stronger than the evidence, and several wordings for one recurring
message. Claims of certainty were checked word by word: every "prove" and
"guarantee" in the interface is a negation ("does not prove translation"),
which is the intended, evidence-limited statement.

**Spelling.** The identifiers are American English (`normalize_sequence`,
`analyze_orfs`, "behavior" in 10 code sites), while the documentation
mixed both systems ("behaviour" 18 times next to "behavior" 11 times,
"labelled", "modelled", "-ise" forms). American English is now used
throughout: 86 words in 8 documents, 9 text sites in 4 source files and
6 in tests, grouped (figures corrected after the closing review; the
first version of this entry said 92 words in 10 documents and 6 code
sites):

| Category | Original | Revised | Where |
|---|---|---|---|
| consistency | behaviour, analysed, labelled, labelling, modelled, normalise(d), normalisation, organises, serialised, recognised, initialised, optimisation, prioritised, parameterised, characterisation, authorisation, ... | behavior, analyzed, labeled, labeling, modeled, normalize(d), normalization, organizes, serialized, recognized, initialized, optimization, prioritized, parameterized, characterization, authorization, ... | all documents in `docs/` |
| consistency | "cut modelled as SpCas9", "Blunt cut modelled", "not modelled", "Modelled cut:" | "modeled" in each; the label keeps its column | `analyses/crispr.py` |
| consistency | "... star activity are not modelled." | "... not modeled." | `core/api.py` (restriction limitation) |
| consistency | "Unrecognised sequence file: ..." | "Unrecognized sequence file: ..." | `fileio/loader.py` |

Words that are correctly spelled with "-ise" in American English
(otherwise, precise, advertise, excise, ...) were excluded by list.

**Recurring messages (8.4).** The interface standard names `Invalid
selection.` as the reply to an unknown menu key; 35 sites used it and two
did not.

| File | Original | Revised | Reason |
|---|---|---|---|
| `analyses/orf.py` (settings menu) | "Invalid option." | "Invalid selection." | one wording per message |
| `analyses/orf.py` (ORF menu) | "Invalid option. Please choose 0-6." | "Invalid selection. Please choose 0-6." | same; the range matches the menu ([1]-[6], [0] Back) |
| `analyses/crispr.py` (guide menu) | "No sequence selected." | "No DNA sequence selected." | the same module says "No DNA sequence selected." elsewhere, and CRISPR targets are DNA |

The primer dialog keeps `ERROR: Invalid selection.`: there the dialog
ends rather than asking again, so the operation did not happen, which is
what `ERROR:` means in the standard.

**Tone and precision.**

| File | Original | Revised | Reason |
|---|---|---|---|
| `main.py` (exit, 2 sites) | "Thank you for using Synapse!" | "Thank you for using Synapse." | professional presentation |
| `analyses/mutation.py` (frame warning) | "frame +1 is probably not the coding frame." | "reading frame +1 is unlikely to be the coding frame." | names the quantity and states the inference as a likelihood (scientific accuracy) |

Left as written after review: "an older RID may simply have expired"
(BLAST) states the ordinary explanation plainly, and the quality bands
"Very Good", "Very Low", "VERY HIGH" are category names, not prose.

Tests: 1 new (`SpellingStandardTests`) checks code and documentation for
31 British forms; it found three sites the first pass had missed (a
docstring in `core/orf_context.py` and a helper name in
`tests/test_stage2_foundation.py`), which were corrected.

### Stage G/H E-values in the BLAST text report (`analyses/blast.py`)

Found while writing the style guide's notation table. The BLAST table and
the HSP lines print E-values with three significant figures (`:.3g`), but
the per-hit block of the text report printed the raw float
(`1.2345678901e-50` for a value that the table shows as `1.23e-50`), and
the evidence note used six (`:g`). All three now go through one helper,
`_evalue_text`, which also prints `n/a` for a missing value instead of
`None`. The configured significance threshold keeps `:g`, because it is a
parameter and should be shown exactly as set.

Tests: 1 new (`test_blast_text_report_shows_e_values_with_three_significant_figures`).
The bundled sample XML has E-value `1e-25`, which both notations print
the same way, so the report part of the test does not distinguish old
from new output; the helper's own assertions (`1.234567e-50` ->
`1.23e-50`, `None` -> `n/a`) do.

### Stage H `docs/STYLE_GUIDE.md`

New permanent reference (Master Prompt section 9): language and spelling,
preferred scientific terms (with "confidence" reserved for statistical
confidence), abbreviations, sequence and variant notation, units,
precision, coordinate conventions, the measured / calculated / predicted /
inferred distinction with Synapse examples, the structure of an
interpretation, message categories and wording patterns for recurring
situations, a summary of the terminal standard, and code documentation
rules.

Every rule in it was taken from the implementation as it now stands, and
writing it surfaced two deviations: the BLAST E-value notation (corrected
above) and the HGVS convention of parentheses for predicted protein
consequences, which Synapse does not follow; the guide records that as
open rather than claiming it.

Its enforcement table names the test for each mechanical rule, and a new
test (`test_style_guide_names_only_tests_that_exist`) fails if the guide
names a test that does not exist, so the table cannot go stale; checked
with a deliberately misspelled test name.

### Stage H README rewritten to section 15.1

The README described running, the API and testing, but not the project's
purpose, its modules, its dependencies, its configuration or its
limitations, and one statement was inaccurate: "BLAST tests that need
BLAST+ or local databases are not part of the default run". They are part
of the default run and are skipped, with the reason shown, when BLAST+ is
absent; the same holds for the Biopython reference tests.

It now covers purpose and status, the eleven menu modules (the planner
marked "not implemented"), installation, a dependency table (runtime:
standard library only; optional BLAST+ and network services; test-only
Biopython and primer3-py; development pyright), terminal and API usage,
the environment variables the code actually reads, the architecture,
testing, the main limitations with figures from `METHODS.md`, and the
documentation set. `LEGACY/` and `Synapse_Website/` are named as not part
of the software. Every statement was checked against the code (the four
CRISPR systems, the default log and export directories, the variables
read through `os.environ`).

### Stage H "Sequence complexity" defined, with its limit (`analyses/sequence.py`, `docs/METHODS.md`)

Audit item P2-5 asked for the "complexity" measure to be renamed. A check
of the four modules that report it (sequence, RNA, CRISPR, mutation)
showed that they compute one and the same quantity: the Shannon entropy
of the A/C/G/T composition divided by 2 bits. That is a composition
measure, and the word "complexity" suggests more: `ACGT` repeated 100
times, a pure tandem repeat, scores the maximum of 1.000.

Rather than rename field names in four modules and the golden master,
the meaning is now stated where readers meet it:

| Where | Original | Revised |
|---|---|---|
| sequence-overview interpretation | "Sequence complexity is high (1.000) based on normalized Shannon entropy of A/T/C/G composition." | "Sequence complexity is high (1.000): the normalized Shannon entropy of A/C/G/T composition. It measures base composition only and does not detect repeats; a tandem repeat of all four bases scores 1.000." |
| `docs/METHODS.md` | no section for the sequence overview | new section: GC/AT content and ranges, skew formulas, the entropy definition and its limit, the Synapse class thresholds, homopolymers, CpG-island criteria (Gardiner-Garden and Frommer 1987), mass conventions, and what is and is not validated |

Category: scientific accuracy. The numbers are unchanged; the
interpretation text is not part of any golden snapshot.

## Phase 6 - Stage I: final verification (2026-11-02)

### External reference comparison: primer Tm against primer3-py (P1-3)

`docs/METHODS.md` and `docs/AUDIT_STATUS.md` still recorded that Primer3
could not be installed, so the audit criterion for P1-3 ("±1 °C against
Biopython/Primer3 at three Mg²⁺ levels") was shown as not run. The Phase 6
initial audit had in fact run 40 comparisons against primer3-py 2.4.0
(max 0.003 °C, `FEATURE_STATUS.md` validation ledger), but neither
document was updated and no test kept the check. (The commit message of
`2851e11` says the criterion "had never been run"; that is corrected
here.) The comparison was now extended and run locally, with no network
access:

- 81 primers: the Primer3 manual example and 18-, 20-, 22- and 25-nt
  windows of a fixed test sequence.
- 36 conditions: Na⁺ 50 and 100 mM; Mg²⁺ 0, 1.5 and 3 mM; dNTP 0, 0.2 and
  0.8 mM; oligo 50 and 250 nM.
- primer3-py with `tm_method="santalucia"`,
  `salt_corrections_method="santalucia"`.
- Result: 2,916 comparisons, max |ΔTm| 0.0032 °C, mean 0.0024 °C. The
  offset is systematic and comes from the gas constant (Synapse 1.98720,
  Primer3 1.987).

The criterion is met with a margin of more than two orders of magnitude.
This shows that Synapse implements the published equations as Primer3
does; it does not show that the nearest-neighbor model predicts
experimental melting temperatures, which no software comparison can.

Test: 1 new (`Primer3ReferenceComparison`, a reduced grid of 31 primers
x 4 conditions, tolerance 0.01 °C), skipped with a reason when primer3-py
is not installed.

### Instability index: the IG dipeptide weight corrected (`analyses/protein.py`)

Found by the new Biopython reference tests (audit Tier 3 item 11, which
asked for the audit's ad-hoc comparisons to become tests that skip
without the reference libraries).

- **Previous behavior.** The dipeptide instability weight for IG
  (isoleucine followed by glycine) was −7.49. Every protein containing IG
  got an instability index lower by 8.49 × 10/L per occurrence. In 31
  random test proteins of 41-400 aa, 20 contained IG and differed from
  ProtParam by up to 1.29 index units.
- **Corrected behavior.** The weight is 1.0, the value in Biopython 1.88
  `Bio.SeqUtils.ProtParamData.DIWV`, which cites Guruprasad et al.
  (1990). (The first version of this entry said the table reproduces
  ExPASy ProtParam; that was not verified.) The rest of Synapse's table equals Biopython's except
  FY (33.6 vs 33.601, a 0.001 difference that is left as is).
- **Reasoning.** The I row read `1.0, 1.0, 44.94, 1.0, -7.49, ...` against
  Biopython's `1.0, 1.0, 44.94, 1.0, 1.0, ...`; -7.49 occurs twice more in
  the same row (IK, IV), which makes a transcription error plausible.
- **Reference.** Guruprasad et al. (1990) through one established
  implementation (Biopython), a single secondary source. The paper's
  table itself was not consulted in this session, because no external
  access was used; checking it is the recommended confirmation.
- **Why the audit missed it.** The two proteins compared with ExPASy in
  the audit (ubiquitin and lysozyme) contain no IG dipeptide.
- **Golden master.** Unchanged: the protein scenario contains no IG.
- **Tests.** `tests/test_external_reference.py` (new, 6 tests): enzyme
  definitions, digests on two topologies and three rotations, the 22 × 22
  end-compatibility matrix, the ProtParam descriptors, the full weight
  table, and a known-answer test for IG that needs no Biopython
  (`AIGA` = 10/4 × (1.0 + 1.0 − 7.49)). On the previous code 20 protein
  subtests, the table test and the known-answer test fail.

### Integration workflows (Master Prompt 10.1-10.4, 10.7, 10.8)

`tests/test_integration_workflows.py` (new, 7 tests) runs complete
workflows through `core.pipeline` and `core.api` instead of single
functions:

- **10.1** A minus-strand ORF flows from detection to selection, RNA and
  protein: its forward-strand coordinates point at the CDS, the RNA is
  the CDS with U, the protein is `MAAAAAAAAAAWKKKKKKKKKK`, and exactly one
  result per module is stored. With NCBI table 2 the internal TGA is read
  as Trp all the way to the protein. Selecting another ORF removes the
  RNA and protein results, the cached sequences and their report entries.
- **10.2** A restriction analysis and a primer design on the same
  sequence share one coordinate frame (EcoRI cut after base 241; every
  primer, reverse primer and amplicon equals its slice of the template),
  and the report carries both with their parameters and the primer
  warnings.
- **10.3** CRISPR guides for a selected ORF lie inside it, match the DNA
  at their 1-based coordinates on either strand, and have their modeled
  cut inside the protospacer; the result records the restriction to the
  ORF and the genome limitation.
- **10.4** A VCF with one matching and one mismatching REF allele: the
  mismatch is reported, both records are kept, and the file's SHA-256 is
  in the parameters.
- **10.7** An exported JSON and Markdown report: the input hash equals the
  SHA-256 of the sequence, every `data_sha256` recomputes from the
  exported data, every entry has parameters, and the Markdown names each
  analysis and the Synapse version.

All seven passed on their first run: they found no new defect, and
the restriction-parameter check passes only since `81ccc55`. The BLAST
(10.5) and literature (10.6) workflows remain covered by the stubbed
client, the offline fixtures and the BLAST+ tests, because a live run
needs external services.

### The CRISPR module could not be imported on Python 3.10 (`analyses/crispr.py`)

`pyproject.toml` requires Python 3.10 or newer, but only Python 3.14 is
installed here, so no test had ever run on 3.10. A static check now
stands in for that run: every source file is parsed with the 3.10
grammar (`ast.parse(feature_version=(3, 10))`), and every f-string is
tokenized to find what PEP 701 first allowed in Python 3.12 (a `#`, a
backslash or the f-string's own quote inside a replacement field).

The check found one site, present since the baseline commit:

| File | Original | Revised |
|---|---|---|
| `analyses/crispr.py:2190` (guide table header) | `print(f"{'#':>3} {'ID':<7} ...")` | `print(f"  # {'ID':<7} ...")` |

Before 3.12, `'#'` inside a replacement field is a syntax error
("f-string expression part cannot include '#'"), so on Python 3.10 and
3.11 the whole CRISPR module, and with it the menu entry, failed to
import. pyright at `pythonVersion = 3.10` does not report this. The
printed header is identical (`'#'` right-aligned in three columns is
`"  #"`).

This is a static check, not a 3.10 run. Running the suite on a real 3.10
interpreter still needs one to be installed, which is the owner's
decision.

Tests: 2 new (`Python310CompatibilityTests`); the f-string test fails on
the previous code and is skipped on interpreters older than 3.12, whose
tokenizer cannot see inside f-strings.

### Documentation corrections after the closing review

`docs/PHASE6_INDEPENDENT_REVIEW.md` (`fc77ff4`) found statements in the
Phase 6 documents that did not match the repository. Corrected, without
any source change:

| Where | Was | Now |
|---|---|---|
| this file (Stage H spelling), `PHASE6_PROGRESS.md` row 35, final report G | 92 words in 10 documents, 6 code sites | 86 words in 8 documents, 9 text sites in 4 source files, 6 in tests (word-level diff of `1a8f2bc`) |
| final report A | Stage I found three defects | items 25-27 come from Stage G (primer coordinates) and Stage I (IG weight, Python 3.10) |
| final report A | 109 tests "at the start" of the previous cycle | at the end |
| final report A | confirmed S items "are corrected" | S9 quantified only, S14 resolved for the API only, S5 addressed |
| `AUDIT_STATUS.md` P1-1 | 312 tests | 314 |
| final report header | attribution "as recorded per commit" | five early commits carry no attribution |
| IG weight (METHODS, FEATURE_STATUS, AUDIT_STATUS, final report, this file) | Biopython's table "reproduces ExPASy ProtParam" | Biopython's table, which cites Guruprasad et al. (1990); one secondary source |
| final report row 26 | "unimportable on 3.10/3.11" | the same, marked as inferred from the grammar, not observed |
| `METHODS.md` references | four entries written without verification | marked as not verified at the source |
| `PHASE6_HANDOVER.md` opening, section 1 | "Phase 6 is not complete" | current status first; the pause text kept as history |
| `PHASE6_PROGRESS.md` "Paused", "Open after the current package", P2-5 row | Stages F-I not started; D17/D18 open; complexity open | marked historical; points to final report J; complexity resolved by documentation |
| final report J | about ten open items missing | added (G·U/RNAfold, primer generations, NEBcutter, real BLAST database, HGVS validator, pipeline printing, complexity names, pyright errors, tool ledger, coverage, two edge-case tests, unverified references, ResourceWarning) |

Not changed, because it is test code and this round changed documents
only: the docstring of
`ProteinReferenceTests.test_instability_weights_match_the_reference_table`
still says "(ExPASy's)", and the module docstring of
`tests/test_external_reference.py` repeats the ExPASy comparison for
ubiquitin from the audit.

## Stage F completion: six modules on the shared renderers (2026-11-02)

The Stage F re-examination (`docs/STAGE_F_UI_REVIEW.md`) found that
Stage F had shared the rules, the `[0]` meanings and the prefix spelling,
but that no module called the renderers in `ui/style.py`. On the owner's
instruction the six modules it compared were converted (`2e47ff9`, `e7eaf58`, `5cb4c9a`, `b181ae3`, `fb4070a`, `fdc968d`, `b80fd37`).

**Method.** Two AST-guided rewriting passes (kept outside the
repository) replaced recurring print patterns at their source positions,
so untouched code, line endings and the two bare-LF lines of `orf.py`
stay byte-identical:

- rule-title-rule triples became `style.header`;
- a heading over a minor rule became `style.section`;
- "Label:   value" prints became `style.field` (a long label or a
  sentence was left alone);
- `print("[n] ...")` runs with their `[0]` entry became `style.menu`;
- `input("Selection: ")` became `ask_choice()`, `input("Label [d]: ")`
  became `ask()`, and the pause prompt became `pause()`;
- prefixed and recognizable messages went through `error`, `warning`,
  `note`, `success` and `invalid_selection`.

What the patterns could not decide was done by hand and is listed in each
commit:

- menu groups for ORF, literature and primer;
- the ORF settings order;
- the ORF candidate table;
- the literature frame;
- the primer Tm dialog;
- the CRISPR title;
- pauses after precondition messages.

**`ui/style.py`** gains `menu()`, `group()`, `field_line()` and
`columns()`, and its docstring now spells "recognizable".

| Before (rendered) | After |
|---|---|
| Literature: unframed banner "Synapse Literature Intelligence v2.1.1" | framed "LITERATURE SEARCH" with module and sources as fields |
| ORF prompt "Select option:", flat menu, no footer rule | "Selection:", groups Discovery / Evaluation / Selection / Settings, shared footer |
| ORF candidate table misaligned, column "Conf." | aligned columns, "Heuristic" |
| ORF settings numbered 1-5, 8, 6, 7 | 1-8 in order (genetic code [6], defaults [7], legacy [8]) |
| CRISPR titled "CRISPR ANALYSIS ENGINE" | "CRISPR GUIDE DESIGN", its main-menu name |
| Context values at columns 9-28 depending on the module | column 24 in all six |
| "SYNAPSE —", "SYNAPSE" title prefixes | none |
| Primer: "ERROR: Invalid selection." ending the Tm dialog | "Invalid selection." and the question again |
| "Run [1] BLAST Search first.", "No variants loaded." scrolled away at once | `NOTE:` and a pause |
| "Export failed: ...", "Exported: ..." without a category | `ERROR:` / `OK:` |

Data, API and golden master are unchanged; no test outside the interface
tests needed a change. Value prompts inside dialogs and the other ten
interface files are listed as open in `TERMINAL_UI_STANDARD.md` section 9.

Tests: 5 new (`RenderedScreenTests`). They render the six first screens
with scripted input and check the frame, title, fields, footer and the
reply to an unknown key, plus a source rule against bypassing
`ui.style`. On the code before the conversion they give 8 failures and
5 errors. They found four more "SYNAPSE" title prefixes the first pass
had missed. The BLAST numbering test now reads the `style.menu()` entries.

## Stage F completion, round 2: the remaining ten interface files (2026-11-02)

The first round (`2e47ff9` ... `623a0a8`) put the six modules named in the
re-examination on the shared renderers. This round does the same for
the other ten interface files (`1ac77a4`, `bf6ae87`, `9b7f2ef`, `85c22ac`, `a4df3d4`, `3853d87`, `69e3f0f`, `7970e5a`, `7ecb88d`). It uses the same rewriting
passes, plus a third one for report line lists. Each commit lists its
visible changes. In summary:

- **Session (`main.py`, `ui/menu.py`).**
  - "MAIN MENU" with labeled groups; "LOADED SEQUENCES" without `[0]`.
  - Loader, type dialog, protein record and error screens without the
    "SYNAPSE -" prefix.
  - The ORF hand-over note no longer claims the RNA codon view uses the
    standard code.
- **Sequence overview, RNA, protein, restriction, cloning, mutation,
  pipeline status.**
  - Titled like their main-menu entries, grouped menus, upper-case
    sections, values at column 24.
  - Context fields on the restriction, cloning and mutation menus.
  - Checkmark lines are `OK:` messages.
- **Tables.** One `columns()` definition each. Four tables were
  misaligned:
  - BLAST database status: the Status header was 42 columns wide, its
    rows 24.
  - CRISPR guide list: the position column was misaligned.
  - Protein composition, both tables: the "3L" column was misaligned.
- **Text reports.** CRISPR, BLAST, primer, variant, literature and
  mutation write their key/value lines through `field_line()`; the CRISPR
  self-test checks the shared form.
- **Prompts and yes/no.**
  - 52 value prompts go through `ask()`.
  - All six yes/no questions go through `confirm()`. This also corrects
    four that the first round had turned into `ask(..., default='y/N')`.
- **Variant menu** in the order of Master Prompt 5.3: normalization and
  filtering before statistics and interpretation. Six entries changed
  number (normalize 8 -> 3, statistics 6 -> 5, analyze 3 -> 6, compare
  5 -> 7, cohort 7 -> 8; filter stays 4).

Not changed, with reasons in `TERMINAL_UI_STANDARD.md` section 9:

- raw line input;
- two CRISPR numeric helpers;
- closing rules;
- analysis result lines;
- the mutation `[0]` deviation;
- the pipeline's printing from `core/`;
- P2-2.

Data, API, golden master and the scientific results are unchanged.
Tests: 3 new and 5 extended in `tests/test_terminal_ui.py`. The rendered
checks cover all eleven menus, the sequence overview and the pipeline
status, and the source rule covers all sixteen files. On the code before
this round they give 11 failures. The extended tests found what the
passes had left: the protein module's own header helper and four
"SYNAPSE ... SELF-TEST" titles, fixed in `7970e5a`. Suite 322/322.

**Correction (closing review, 2026-11-03).** The summary above
overstated round 2 in two points:

- **"Data, API, golden master and the scientific results are
  unchanged" is true, but the workflow was not.** From round 1 until
  `1516339`, Enter on fifteen bracketed prompts returned the hint text
  (next entry).
- **"All sixteen files on the shared renderers" is true for adoption,
  not for completeness.** 59 direct prints remained, none of them
  documented:
  - 31 messages without a category;
  - 19 headings printed directly;
  - 9 dynamic option lists in BLAST and CRISPR.

  The rendered tests drive the first screen of each menu: 17 of the
  38 `style.menu` call sites.

`docs/STAGE_F_FINAL_REVIEW.md` has the measurements.

### Fix: pressing Enter in six dialogs broke them (Stage F regression, 2026-11-03)

Reported by the owner: a literature search never ran. The Stage F
rewriting pass had turned prompts of the form `Label [hint]: ` into
`style.ask("Label", default="hint")`, so Enter returned the hint text
instead of an empty answer. Where the hint was the value itself
(`[10]`, `[30]`, `[75]`) that was harmless; where it was a description,
the dialog broke:

| Dialog | Enter returned | Effect |
|---|---|---|
| Literature: From / To date `[none]` | `none` | "Invalid configuration: date_from must use YYYY ..."; no search ever ran |
| Primer: Na+ `[50.0 mM]`, Mg2+ `[1.5 mM]`, primer `[250 nM]` | the hint | `float()` failed; primer design stopped at the conditions |
| BLAST (local): word size `[default]`, filtering `[default/none]` | the hint | an error instead of the default settings |
| Cloning planner: strategies `[default 5]` | `default 5` | an endless "Please enter a positive integer." |

Introduced in `5cb4c9a` (BLAST), `b181ae3` (literature), `fb4070a`
(primer) and `85c22ac` (cloning). All fifteen converted prompts keep
the bracket in the label again, so Enter returns `""` and the caller
applies its default, as before Stage F; the prompt text is unchanged.
`main.py`'s report path keeps a real default, which the caller
already used for an empty answer. The literature configuration error
is now an ERROR message.

No test had driven these dialogs with an empty answer. New
`EmptyAnswerTests` do: the literature dialog (dates None, 30 results,
relevance), the primer conditions (50, 1.5, 0.2, 250), the cloning
strategy count (5), and a source rule against `ask(..., default=...)`
in the modules. All four fail on the previous code.

## Stage F closing review, recommendations 1-3 (2026-11-03)

`docs/STAGE_F_FINAL_REVIEW.md` measured what the completion rounds had
left; the owner released recommendations 1-3.

- **Claims (`46710cd`).** The final report, the UI review and the
  progress log no longer call Stage F "unified" for round 2. They now
  separate three things: the partial standardization of Stage F, the
  adoption reached in the completion rounds, and what was still open.
  They also state that the workflow was not preserved until `1516339`.
- **Messages (`4e2cace`).** 44 plain prints go through the helpers of
  section 5 of the standard:
  - an invalid choice gets `invalid_selection()`;
  - a failed input check is an ERROR;
  - a cancelled action or context is a NOTE;
  - a completed action is OK;
  - BLAST's "IMPORTANT:" is a WARNING.

  Control flow is unchanged.
- **Headings (`e559255`).** The primer report sections, the literature
  WARNINGS / SOURCE ERRORS / ABSTRACT / KEYWORDS, the CRISPR ranking and
  the BLAST privacy and resume blocks go through `section()`. Each now
  carries a rule under the title; blank lines are kept. "BLAST SEARCH
  COMPLETED" is an OK message.
- **Option lists (`b9ec344`).** Six BLAST lists and three CRISPR lists
  go through `style.menu(zero="Cancel")`. They gain a group label, the
  blank line before `[0]` and the closing rule. The hit and record lists
  show `[0] Cancel` instead of "(0 cancel)" in the prompt. The BLAST
  evidence profile, a result list, is numbered "1.".
- **Also corrected:**
  - three group labels with a stray `"\r\n"`;
  - the pointer to "[3] BLAST Database Manager", which is [4];
  - the CRISPR batch entry label.
- **Tests (`73a1d93`, `ccc987b`).** `DirectOutputTests`; 68 sites on
  `1516339`, none now.
- **`FEATURE_STATUS.md` (`163c54d`).** Variant, ORF-settings and BLAST
  menu numbers as in the code; the summary gives 328 tests and the final
  defect status.

No scientific calculation, data model or golden-master value changed.
Suite 328/328, golden master 12/12, pyright 112 with no new diagnostic.

## Stage F closing review, recommendations 4 and 5 (2026-11-03)

- **core/ui (`b25e800`, `b4e5a53`).**
  - `PipelineContext.status_report()` returns the status as data.
  - `print_status(render)` takes a `StatusRenderer`; `[P]` passes
    `ui.report.print_pipeline_status`.
  - core imports `ui` nowhere; the screen is byte-identical.
  - `LayeringTests` runs core in an interpreter that refuses `ui`.
  - API change: `print_status()` now needs the renderer argument. Its
    only callers were `main.py`, `SynapsePipeline` and a test.
- **Single-Enzyme Analysis (`8eddc72`).**
  - Previous behavior: every enzyme raised `KeyError: 'left'`, since at
    least the audit commit, because the result went to the two-enzyme
    report.
  - Corrected behavior: a report of its own from the fields
    `analyze_cloning()` returns. No calculation changed.
- **Tm model dialog (`c78800e`).**
  - Previous behavior (since `fb4070a`): Enter answered "Invalid
    selection.".
  - Corrected behavior: Enter chooses nearest-neighbor, as the prompt
    says.
- **Messages (`61d0d38`).** Two duplicate-enzyme errors in cloning and the
  CRISPR export-extension error are ERROR lines.
- **Rendered menus (`a39c39f`).** `tests/test_rendered_menus.py` drives all
  47 menus: entries, numbering, `0`, invalid input, error states and
  result views. The NCBI client is mocked.

Golden master unchanged; suite 356/356; pyright 111 with no new
diagnostic.

## Phase 6 completion run, packages A and B (2026-11-11)

Each commit message holds the measurements and the tests that fail on the
previous code.

- **HTTP error responses (`68e59d1`).** The literature and BLAST clients
  left the response an `HTTPError` carries open (the ResourceWarning in a
  retry test); they close it now.
- **Mutation prompt and messages (`8379ac2`).** The sequence prompt has
  `[0] Cancel` (`SequencePromptCancelled`); ten messages that were drawn as
  fields are ERROR/NOTE lines, with a rule that finds them; `MissingFileTests`
  (Master Prompt 11.4).
- **Cloning planner (`b4b16bb`, P2-6).** Each enzyme analysis is computed
  once per search: 21.7 s -> 3.2 s, results byte-identical.
- **Type IIS (`4939688`, P2-7).** BsaI, BsmBI, BbsI and SapI are searched on
  both strands and cut outside their sites; 140 cuts against Biopython 1.88,
  0 differences. Not used in cloning; Golden Gate planning not implemented.
- **HGVS-style notation (`b7eb1de`).** Predicted protein changes are written
  `p.(R12C)`, `p.(R12=)`, `p.(W24*)`, `p.(M1?)`, `p.(*110Qext*?)`; used in the
  mutation report and in the variant field `protein_change` (it held the
  consequence class). Not checked with an HGVS validator.
- **Restriction coordinates (`1ae2833`).** Unsuffixed `site_position`,
  `top_cut`, `bottom_cut` are 1-based in every record; `*_0based` keys give
  0-based values (maps used 1-based, digests 0-based before).
- **RNA hairpins (`fedb750`, P2-5).** G·U pairs are counted as a separate
  figure (golden reference: 58 Watson-Crick, 260 with G·U).
- **ORFs (`f020a1d`, P1-9).** ORFs across the origin of circular sequences
  and optional 5'-partial ORFs; rotation-invariance test.
- **Primer structure screens (`27255df`).** Compared with primer3-py ΔG:
  Spearman ρ −0.78 / −0.74 / −0.51; classes overlap (METHODS).
- **Variant cross-check (`2c68e61`, P2-4).** `core.api.cross_check_variants`.
- **Alignment (`f9d4356`, `e5bb456`, `15e2ff1`, P2-6).** One traceback byte
  per cell (46 -> 1.1 bytes); exact banded alignment of similar sequences
  (10 kb 23 s -> 0.18 s); banded Levenshtein similarity (10 kb 0.10 s);
  faster hairpin screen (primer design 10 kb 135.6 s -> 25.6 s). Results
  identical.
- **Records (`b8e92d8`, `690c9e8`, `bd91a24`).** Continuation record and the
  complete handover.

## Phase 6 completion run, second part (2026-11-12 and 2026-11-13)

- **Approval policy (`bf4d096`).** The owner adopted the Master Prompt's
  approval rules (1.3, 14.5). The handover quoted a commit subject with
  British spelling, which made the suite fail at `bd91a24` (420/421); fixed.
- **Primer design 10 kb (`aba3fb0`, P2-6).** 71 s -> 10 s CPU time on the
  same machine under the same load, results identical: a Tm-tier screen
  with an exact upper bound of the screening score decides which windows
  are measured in full; the dimer stem is a longest common substring; the
  hairpin, repeat and 3' screens use string operations; pure functions are
  cached. Evidence: `tests/test_primer_screen_equivalence.py` (the replaced
  algorithms as references, screen vs. full measurement, the bound).
- **D20, D21 (`c6f0df9`).**
  - D20, previous: `_v5_candidate_quality` read hairpin and self-dimer
    stems from fields the thermodynamic profile does not have, so its
    documented structural penalty was always 0. Corrected: the stems come
    from the candidate. Some designs change (third pair of a 900-bp
    template); golden master unchanged.
  - D21, previous: `tm_method="nearest-neighbor"` passed validation and
    then raised "Unknown Tm method". Corrected: one alias set.
- **Primers around a restriction site (`99150da`, Master Prompt 10.2).**
  `restriction_site_choices`, `design_primers_for_restriction_site`,
  `core.api.design_primers_for_restriction_site`, primer menu entry 2
  (Single Primer Analysis moved to 3), pipeline parameters. Each pair
  carries the predicted digest of its amplicon.
- **P2-2 (`50dea66`).** The dialogs of eleven analysis modules (130
  definitions, 7 aliases) moved byte for byte to `ui/<module>_menu.py`;
  `core/textformat.py` holds the string builders, re-exported by
  `ui.style`; `SynapsePipeline.select_orf(chooser=...)` replaces
  `interactive=True` (no caller). API change: the dialog functions are
  imported from `ui.<module>_menu`, no longer from `analyses.<module>`.
  Evidence: 55 rendered screens byte-identical; no analysis module imports
  `ui`, checked in a subprocess with `ui` unimportable.
- **Pipeline parameters (`b1a012f`, T4).** The sequence overview and the
  three cloning dialogs record their parameters.
- **D22, D23, typing (`757958f`).**
  - D22, previous: `mutation.orf_reference_cds` raised AttributeError for
    an ORF without a sequence. Corrected: MutationAnalysisError.
  - D23, previous: the ranking's risk flags, structural warnings,
    screening scores and model Tm difference were attributes outside the
    dataclass fields and never reached JSON or Markdown exports. Corrected:
    dataclass fields with the defaults the ranking assumed.
  - `calculate_tm_nearest_neighbor` overloads; pyright 108 -> 58.
- **Coverage, D24 (`ae87adb`).** coverage.py measured 68 % -> 70.9 %
  (`docs/TEST_COVERAGE.md`); the protein and RNA views and eleven
  never-called functions are executed with output checks.
  D24, previous: the BLAST text export and the local query FASTA were
  written with CR CR LF line ends. Corrected: `newline=""`.
- **Ledger, Python 3.10, pUC19 (`bf4b9ee`).** `docs/TOOL_AND_DATA_LEDGER.md`;
  the suite runs on Python 3.10.11 (project-local embeddable package);
  pUC19 L09137.2 is a fixture and its 132-digest comparison with Biopython
  runs in the suite.
- **BLAST dialogs (`dabc2f1`, T11).** Previous: the hit selector answered
  2 of 1 with "list index out of range" and showed a hit from the end of
  the list for −1; the query-type question cancelled on an unknown key.
  Corrected: only listed numbers; unknown keys are answered and re-asked.
  The literature limitation line is a NOTE.
- **README and records (`a65cb04`).**
- **D25 (`a5e5a94`).** Previous: one network error while polling a
  submitted NCBI search ended it, and the message did not name the RID.
  Corrected: up to 3 consecutive transient errors are retried; a failed
  submitted search names its RID. Found in a live search; the repeated
  live search (pUC19 401-700 against `nt`) succeeded. The primer3-py
  comparisons are skipped with the reason when Windows application
  control blocks its compiled module (observed on 2026-11-13, temporarily).
- **Variant check in the menus (`cdbf69c`, T12, P2-4).** The variant
  dialog keeps the variants of the contig mapped onto the loaded sequence;
  the primer and CRISPR designs show "Variants in binding sites".

Golden master unchanged throughout; suite 463 tests, all passing (the
opt-in benchmark skipped); pyright 58 errors, nothing new against the
`ab4bd2a` baseline.

## Continuation run (Master Execution Prompt, from 2026-11-13)

Governed by the Phase 6 Master Prompt and the owner's Master Execution
Prompt of 2026-11-13. History and binding decisions:
`docs/PROJECT_HISTORY_AND_DECISIONS.md`.

### AI Experiment Planner: COMING SOON (`main.py`, `ui/menu.py`)

- Previous: `[11] AI Experiment Planner (not implemented)`; selecting it
  printed "not implemented ... deferred" (D15, 2026-10-31).
- Now: `[11] AI Experiment Planner (Coming Soon)`; selecting it prints
  "AI Experiment Planner: COMING SOON. This feature is planned but not yet
  available. Nothing is run or stored." Nothing is implemented behind the
  entry; the numbering is unchanged.
- Reason: owner decision of 2026-11-13 (Master Execution Prompt §3). The
  wording keeps D15's requirement that the entry must not suggest a
  working feature.
- Tests: `PlannerPlaceholderTests` drives `main.py` and checks the entry.

### Optional energy-based RNA structure with ViennaRNA (P2-5)

- New: `analyses/rna_structure.py`, `core.api.rna_secondary_structure`,
  RNA menu entry `[6] Secondary Structure (ViennaRNA, optional)`. With the
  ViennaRNA package (2.8.0, installed project-locally as a wheel; no
  system install) Synapse computes the MFE structure (Turner 2004, 37 °C)
  and annotates each hairpin-screen candidate with its free energy alone
  and its support in the MFE structure. Without the package the view and
  the API say "not available" with the reason.
- Previous: RNAfold was recorded as BLOCKED (assumed to need a
  system-wide binary).
- Comparison: the hairpin screen finds 64 of 74 MFE hairpins (86 %), but
  only 35 % of its 3,307 candidates are favorable on their own
  (`docs/METHODS.md`). No existing result changes: the screen, the RNA
  analysis output and the golden master are untouched.
- Tests: `tests/test_rna_secondary_structure.py` (12; the ViennaRNA tests
  skip without the package).

### Variant normalization compared with `bcftools norm`; reference-only records reported

- Comparison (no behavior change): Synapse's split-and-normalize result
  equals `bcftools norm -m-` on the bcftools 1.25 test vectors for 29 of
  29 alleles; the REF mismatch that `-cx` excludes is the one Synapse
  flags. Fixtures with SHA-256 and licence in
  `tests/fixtures/bcftools_norm/`.
- Previous: a VCF record without an ALT allele (ALT `.`) produced no
  variant and was dropped without notice.
- Now: it still produces no variant (it is a reference-only site), but
  `parse_vcf(..., skipped=[])` lists it, `core.api.analyze_vcf` warns and
  returns `data["reference_only_records"]`, and the variant menu prints a
  NOTE after loading (Master Prompt 10.8, no silent loss).
- Tests: `tests/test_bcftools_norm_reference.py` (7).

### ORF search compared with NCBI ORFfinder

- Comparison (no behavior change): ORFfinder's web service on pUC19
  (L09137.2), genetic code 11, four settings; all 69 ORFfinder ORFs are
  reported by Synapse with identical strand, start, stop, length and
  partial status. The result tables are stored as a fixture.
- Correction of a documentation claim: METHODS said the 75 nt default
  "follows ORFfinder". ORFfinder's minimal length excludes the stop codon
  and Synapse's includes it, so Synapse also reports ORFs of exactly 75 nt;
  `min_length=78` reproduces ORFfinder's 75. The default is unchanged.
- Tests: `tests/test_orffinder_reference.py` (4).

### HGVS DNA descriptions in the mutation report (checked with Mutalyzer)

- Previous: `mutation_hgvs_c` prefixed Synapse's own notation with "c.",
  and the report printed that notation: an insertion before base 11 was
  `11insTT` (HGVS needs both flanking positions, and here the change is a
  duplication), deletions repeated the deleted bases, replacements were
  `20_21GT>CC`, and no 3′ rule was applied.
- Now: `mutation_hgvs_c(mutation, reference_sequence)` writes HGVS syntax
  (`c.7C>T`, `c.16del`, `c.10_11insTT`, `c.9_10dup`, `c.11_12delinsCCC`)
  and, with the reference, applies the 3′ rule and recognizes
  duplications; `hgvs_c_allele` describes all changes as one allele and
  merges changes that touch into one delins. The report shows both.
- Evidence: Mutalyzer 3 `description_extract` on 11 synthetic sequence
  pairs, 11/11 identical (one in the equivalent dup form); the p. forms
  parse in Mutalyzer. Mutation positions, types, effects and the golden
  master are unchanged; only the printed notation changed.
- Tests: `tests/test_hgvs_mutalyzer_reference.py` (5).
- Record: the test totals in the commit messages of `70ff2b0` (484) and
  `a1e3520` (488) are one too high; the suites had 483 and 487 tests.

### Instability-index weights compared with ExPASy ProtParam

- Comparison (no behavior change): all 400 dipeptide weights were
  recovered from 420 ProtParam measurements on synthetic probes and
  compared with Synapse's table. 399 agree; KQ differs (ProtParam 24.68,
  Synapse and Biopython 24.64). The table is unchanged because the 1990
  paper is needed to decide which value is the published one.
- Tests: `tests/test_instability_expasy_reference.py` (3), with the
  measurements as a fixture.

### NCBI `core_nt` offered for remote BLAST (2026-11-17)

- Previous: the remote nucleotide databases were `nt`, `refseq_rna` and
  `refseq_genomic`; NCBI's default since 2024, `core_nt`, was missing.
- Now: `core_nt` is listed first among the remote nucleotide databases
  and is in the `update_blastdb.pl` catalog (not recommended for
  download). `nt` and the default of the API (`database="nt"`) are
  unchanged, so scripts behave as before.
- Menu: in the remote nucleotide list, `core_nt` is entry 1 and `nt`
  moves from 1 to 2.
- Evidence: live search against `core_nt` (RID FR8EK6FN014, public pUC19
  fragment): NCBI reported the database as `core_nt` and gave the same
  top hits as `nt`. Tests: `tests/test_blast_core_nt.py` (4); the two
  rendered-menu assertions that list the databases were updated.

### Dialog tests and message prefixes (2026-11-17)

- New `tests/test_dialog_paths.py` (17): the variant, mutation, ORF and
  literature dialogs driven through missing and invalid input, cancel,
  settings and exports (`docs/TEST_COVERAGE.md`: 23-38 % -> 72-87 %).
- Message prefixes (Master Prompt 5, 8.4): "Run a sequence comparison
  first." (mutation, twice) is now a NOTE; the coding-frame warning is one
  `WARNING: Coding frame assumption: ...` line per finding instead of a
  hand-written "WARNING - coding frame assumption:" list; "Search did not
  produce a run.", "Enter a valid result number." and "Article number is
  outside the result range." (literature) are ERROR lines.
- ORF settings: a second note says that Synapse counts the stop codon in
  the minimum length and ORFfinder does not (78 nt here = ORFfinder's 75).
- No analysis result changed; golden master unchanged.

### Documentation reconciled with the continuation run (2026-11-17)

- README: the ORFfinder limitation now states the comparison and its
  scope; the hairpin screen and HGVS numbering limits are named; the test
  count is current.
- FEATURE_STATUS: test count 516; the instability index records the
  ExPASy comparison; the primer structure screens record the primer3-py
  comparison of 2026-11-11 (the row still said "Not validated"); PubMed,
  Europe PMC and OpenAlex record the live check of 2026-11-12 (Crossref
  was not part of it); origin-spanning ORFs are named as outside the
  ORFfinder comparison.
- AUDIT_STATUS: P1-1 test counts, P1-9 (ORFfinder), P1-10 (HGVS) and
  open validation 9 (ExPASy, KQ open) updated.
- CFD scoring (FEATURE_STATUS, AUDIT_STATUS) is marked deferred under the
  owner's decision of 2026-11-17; its content is unchanged.
- PROJECT_HISTORY_AND_DECISIONS: that decision, and the status of every
  previously open item after this run.

## Phase 7 – engineering stabilization (from 2026-11-17)

Governed by the owner's Phase 7 Master Execution Prompt (engineering
quality only; no new scientific features; CRISPR/CFD excluded; the AI
Experiment Planner stays Coming Soon).

### Package A: FEATURE_STATUS checked row by row

- Method: every row's named functions were looked up in the code, and
  their bodies were checked against line coverage of the full suite.
- Out of date and corrected: the RNA codon view (uses the ORF's genetic
  code since P1-9, not "always table 1"); the cloning digest fragments
  (S17 corrected long ago); local BLAST (blastp works since D14-F2); the
  core API (all entry points present since S14); custom BLAST databases
  (tested with real `makeblastdb`); seven views marked "Reachable" that
  tests drive; the literature menu entry (live-checked); the size-guard
  note (10 kb similar sequences now fast).
- Overclaimed and now tested: the double digest,
  `diagnose_thermodynamic_failure`, `normalize_indel_left`,
  `compare_proteins` (no menu caller) and the Europe PMC, OpenAlex and
  Crossref providers. FEATURE_STATUS said "Execution-tested" but no test
  ran them. New `tests/test_library_functions.py` (14).
- Deleted names removed from the table: `_effective_monovalent_salt_mM`
  and `_professional_failure_diagnostics` (D16 part 2).
- Section headers no longer carry line counts, which went stale with
  every change. The summary and the validation ledger (§17: 24
  validations, current open list) are brought up to date.
- New `tests/test_documentation_consistency.py` (2): every function the
  inventory names must exist in the code.
- Addenda added to `PHASE6_FINAL_REPORT.md` and
  `PHASE6_FINAL_ASSESSMENT.md`; their original text is kept.
- Code: unchanged. A dead `elif pmid` branch in the Europe PMC URL
  construction was noted, not changed.

### Package B: pyright 58 -> 11 (the 11 are in the excluded CRISPR module)

Each finding was traced to its cause; no check was disabled, no broad
`Any` was added, and one targeted ignore carries its reason.

- `analyses/literature.py`: `relevance_components` was declared
  `dict[str, float]` but also holds term lists and the ranking-method
  text; it is now a `TypedDict` (`RelevanceComponents`) with the real
  keys. The date-range check, the Crossref DOI quote and the header-based
  query suggestion use local values pyright can see are set; the self-test
  compares scores defensively. Same behavior.
- `analyses/orf.py`, `analyses/rna.py`: the Shine-Dalgarno search sets
  motif and position together, and the early return now checks both;
  `get_rna_complement` is typed `str -> str`; the ORF nesting record uses
  one `nesting` tuple instead of two variables only bound on one branch.
- `analyses/protein.py`: the mass-range bounds are checked together.
- `analyses/blast.py`: the runtime alias directory is `None` unless an
  alias was built, and the environment uses it instead of re-evaluating
  the condition.
- `core/report.py`: `to_dict` is looked up once with `getattr`.
- Tests: `assert x is not None` before using results that can be `None`;
  real `email.message.Message` headers for `HTTPError`; a package import
  of the genetic-code fixture instead of a `sys.path` change; a typed
  cast for the primer stand-ins; one `pyright: ignore[reportArgumentType]`
  where Biopython's stub declares `table: str` although its
  documentation and code accept NCBI table ids.
- Left: 11 findings in `analyses/crispr.py`, not edited because CRISPR
  work is excluded in this phase.
- Suite, golden master and Ruff (pyflakes) unchanged.

### Package D, part 1: BLAST dialog tests; Remove Database left a file behind

- Defect found by the new test: "Remove Database" in the BLAST database
  manager left the `.njs` file (the JSON metadata of BLAST+ version 5
  databases, which `makeblastdb` writes) because the manager's list of
  database files did not contain `.njs`/`.pjs`. The same list feeds the
  file count and size in "Database Information", which were too low.
  Corrected in `LocalBlastDatabaseManager._expected_extensions`; the test
  failed before the change and passes after it.
- New `tests/test_blast_dialog.py` (6): every result view after two
  stubbed searches (configuration, summary, hit detail, evidence,
  taxonomy, comparison, report, self-test, catalog, status), the JSON,
  CSV and text exports, missing and unsupported export paths, a declined
  remote search, the RID resume (missing, declined, resumed, failed) and
  the database manager on a real `makeblastdb` database (create, inspect,
  verify, wrong confirmation, remove; skipped without BLAST+).

### Package D, part 2: the session, the single-primer analysis and the long-sequence overview

- New `tests/test_main_session.py` (7): `main.main()` driven in-process
  (the subprocess smoke test is not measured by coverage). It covers a
  missing file, the full DNA session with reports in JSON and Markdown,
  protein before an ORF, the planner entry (Coming Soon), every module
  entry except CRISPR, sequence and file switching, a protein record, an
  undetermined type declared or cancelled, an invalid record, and
  `report_error`. `main.py` coverage 34 % -> 83 %.
- New `tests/test_primer_and_overview_dialogs.py` (4): the single-primer
  analysis (empty, invalid, cancelled, bad conditions, nearest-neighbor
  and Wallace) and the overview profiles of a 2.7 kb sequence.
- Found and left for package C: `print_protein_report` (about 600 lines)
  and `print_primer_report` have no caller.

### Package D, part 3: literature search pipeline; BibTeX export escaped too little

- Defect found by the new test: the BibTeX export escaped only `\`, `{`
  and `}`. A title such as "Lactase & diet" or "50% of cases" produced a
  `.bib` entry on which LaTeX stops ("Misplaced alignment tab character").
  `&`, `%`, `$`, `#`, `_`, `~` and `^` are now escaped in title, authors,
  journal, abstract and note; the DOI field keeps its characters (it is an
  identifier). The tests fail on the previous code.
- New `tests/test_literature_pipeline.py` (8), HTTP replaced by
  recorded-shape stand-ins: PubMed search (structured abstract,
  collective author, DOI and PMC ids, date filters, no hits, malformed
  XML, invalid arguments), the combined search (one provider failing,
  open-access and date filters, DOI de-duplication, canonical links),
  metadata enrichment from Europe PMC, PubMed and Crossref (and failures
  kept as warnings), and the five export formats.

### Package C, part 1: Ruff default rules clean; Ruff configured

- Ruff's default rule set (pyflakes plus pycodestyle E4/E7/E9) reported
  20 findings. Fixed 19 without behavior change: 15 one-line `if`/`for`
  statements in the literature exports and summary split onto two lines;
  a second, redundant `from dataclasses import ...` in the middle of
  `analyses/primer.py` removed (the names were already imported at the
  top); the ambiguous names `l` renamed (`leu` in the aliphatic index,
  `text` in the GenBank reader and one test).
- `pyproject.toml` now has a `[tool.ruff]` section (target Python 3.10,
  legacy folders excluded). The remaining finding, a lambda assignment in
  a CRISPR test, is ignored for that file with the reason stated, because
  CRISPR code and tests are not changed in this phase.

### Package C, part 2: `main()` from 632 to 383 lines

- The sequence overview computed 16 values into separate variables and
  passed them to `print_results` in two identical 20-argument calls (on
  selection and for menu entry 1), and once more into the stored result.
  `compute_sequence_overview(dna)` now returns them by name and
  `show_sequence_overview()` draws them; the stored result is the same
  dictionary.
- The ORF step (analysis, pipeline storage, confirmation screen), the
  protein step from the selected ORF, its "no ORF selected" message and
  the AI Experiment Planner message (Coming Soon, unchanged) are separate
  functions. The CRISPR branch is unchanged.
- Evidence: two scripted sessions (the example file through every module
  entry except CRISPR, and pUC19 through ORF, protein, RNA and status)
  give byte-identical output before and after (64,706 and 26,122
  characters); `test_main_session`, the smoke tests and the golden master
  pass.

### Packages E and F: continuous integration, dependencies, installation

- New `.github/workflows/ci.yml` (GitHub Actions, no paid service):
  standard-library install and test suite on Python 3.10-3.14 (Linux)
  and 3.14 (Windows) with the golden master as its own step and the
  skipped tests listed with their reasons; a reference job with Biopython,
  primer3-py, ViennaRNA and BLAST+ in which only the opt-in benchmark may
  skip, with a coverage report; a lint job with Ruff and the pyright gate.
  The repository has no remote, so the workflow has not run on GitHub;
  every command it runs was run locally (YAML parsed, scripts executed,
  clean install below).
- New `scripts/ci/pyright_gate.py`: fails on any pyright error outside the
  documented allow-list (11 findings in `analyses/crispr.py`, excluded
  from changes) and asks for the number to be lowered when it drops.
- New `scripts/ci/check_skips.py`: lists skipped tests with reasons (also
  in the GitHub job summary); with `--only` it fails on any other reason.
- `pyproject.toml`: the `dev` extra lists the tools the project uses
  (pyright 1.1.414, Ruff 0.14.4, coverage >= 7.10) instead of pytest,
  which no test imports and no instruction used; the unused
  `[tool.pytest.ini_options]` is removed (the suite stays runnable with
  pytest). Python version and topic classifiers added. Runtime
  dependencies remain none.
- Clean-environment check: a new virtual environment with nothing but
  `pip install .` installs Synapse without dependencies; the modules
  import from `site-packages`, the golden master passes and the suite
  passes with 22 visible skips (Biopython 11, ViennaRNA 6, primer3-py 4,
  benchmark 1).
- `.gitignore`: build artifacts (`/build/`, `/dist/`, `*.egg-info/`),
  `.coverage` and `.ruff_cache/`.
- README: development setup in a project-local environment and the
  checks the CI runs.

### Package G: messages and prompts

- Five hand-padded lines in the single-primer analysis ("ΔH:", "ΔS:",
  "ΔG at 25 °C:", "3' GC clamp:" twice) are drawn with `style.field()`, so
  their values align with the other fields. TERMINAL_UI_STANDARD had said
  no such line remained; corrected there.
- BLAST: an unsupported export suffix is an ERROR line; "No hits were
  returned" is a NOTE; three database-manager prompts no longer begin
  with a carriage return (a blank line is printed before them instead).
- Restriction/cloning: "Cloning module could not be imported." is an
  ERROR line. Literature: the date-validation remark is a WARNING line.
- Main session: "Protein Analysis will not use the full DNA sequence as a
  fallback." is a NOTE like the two lines above it.
- New `tests/test_message_rules.py` (3): no hand-padded "Label:   value"
  print and no prompt beginning with a line break in the dialogs (CRISPR
  excluded); both rules fail on the previous code.
- BibTeX export escaping (earlier in Phase 7, `ca5d53a`) belongs to this
  package as well.

### Package H: the genetic code in the RNA and protein report entries

- Previous: the stored RNA and protein results recorded only their source
  as parameters; the genetic code they used was inside their data, not in
  the parameters a reader of the report looks at.
- Now: `PipelineContext.store_rna_analysis` adds `input_mode` and
  `genetic_code_table`, `store_protein_analysis` adds `genetic_code_table`
  and `reading_frame`, each read from the result and only when present
  (a protein record has no translation and keeps only its source). No
  calculation changes.
- New `tests/test_report_traceability.py` (5): a table-11 session on pUC19
  exported to JSON; every entry carries versions, time, hashes,
  parameters and limitations, and ORF, RNA and protein all state table
  11. It fails on the previous code.
- METHODS, "Reproducibility": what the report records, per header and
  entry, and where external tools record their versions.

### Package I: licence recommendation and governance

- New `docs/LICENSING_AND_GOVERNANCE.md`: why a licence is needed, the
  licences of every third-party item in or around the repository,
  permissive, copyleft and source-available options with their
  consequences, and the owner decisions they require (licence, copyright
  holder, citation file, two unused legacy BLAST result files, two
  unreachable report printers). No licence was selected and no LICENSE
  file was added.
- New `CONTRIBUTING.md`: setup, the checks to run, and the project rules
  (golden master, no fabricated data, labeled heuristics, layers,
  language, optional dependencies, documentation, history, owner
  decisions).
- README: a short licence section that points to both.
- Repository hygiene checked: no generated file is tracked.

### Package J: final verification

- Suite 565 tests, 0 failures on CPython 3.14.6 (5 skipped) and 3.10.11
  (23 skipped, reference packages not installed for 3.10); golden master
  12/12; pyright gate passes (11 errors, all in `analyses/crispr.py`);
  Ruff clean; coverage 81 %; clean-environment install and wheel checked;
  `git diff --check` clean; ReadOnly 15/15.
- FEATURE_STATUS and AUDIT_STATUS carry the final test count; FEATURE_STATUS
  records the BibTeX escaping and the new "Remove database" row.
- Reports: `docs/PHASE7_FINAL_REPORT.md` (A-L), `docs/PHASE7_CONTINUATION.md`.

## Phase 7.1 – architectural modularization (from 2026-11-18)

Governed by the owner's Phase 7.1 prompt (internal architecture only; no
new scientific features; CRISPR/CFD excluded; the AI Experiment Planner
stays Coming Soon). Audit: `docs/PHASE71_ARCHITECTURE_AUDIT.md`.

### Package A: baseline and characterization tests (`cd23f74`)

- Baseline at `9077f1e` verified again (565 tests; primer3-py now loads on
  this machine, so its four comparisons ran and passed).
- `tests/test_module_namespaces.py`: the 264 top-level names of
  `analyses.blast` and `analyses.primer` stay reachable.
- `tests/test_characterization.py` with
  `tests/fixtures/characterization/blast_primer.json`: 78 scenarios
  recorded from `9077f1e` and compared exactly (floats by `repr`).
- Code: unchanged.

### Package B: `analyses/blast.py` split into `analyses/_blast/`

- **Before.** One module of 4,667 lines with ten responsibilities.
- **After.** `analyses/blast.py` (217 lines) is the public module. It keeps
  its docstring and the `python -m` entry point and re-exports all 125
  names it defined. The code is in 14 modules of `analyses/_blast/`:
  `paths`, `config`, `models`, `query`, `parsing`, `ncbi`, `executables`,
  `runtime`, `databases`, `local_engine`, `evidence`, `reporting`,
  `workflow`, `selftest`. The import graph between them has no cycle.
- **How.** A script cut the module at statement boundaries and generated
  the imports from the names each piece uses. An independent check found
  all 127 statements verbatim in their new place.
- **The only edit inside moved code.**
  - Four `Path(__file__).resolve().parent[.parent]` lookups (database
    root, executable search, user paths, one self-test) now read
    `ANALYSES_DIR[.parent]`.
  - `ANALYSES_DIR` (`analyses/_blast/paths.py`) is the `analyses`
    directory, which is what `Path(__file__)` meant in the old location.
  - The snapshot records the database root, the project-local executable
    directories and the user-path base: all unchanged.
- **`rank_hits`** moved to `models` next to `BlastRun`, because
  `BlastRun.ranked_hits` calls it. Placing it with the evidence code would
  have created an import cycle.
- **Names no longer in the namespace.** The standard-library and
  `core.alphabets` / `core.textformat` names that `analyses.blast` only
  imported (`json`, `time`, `urlopen`, `Path`, `RNA_IUPAC`, `minor_rule`,
  ...) are no longer attributes of it.
  - Two tests patched them: `B.time.sleep` in `test_blast_polling` and
    `B.urlopen` in `test_stage3_additions`.
  - They now patch `analyses._blast.ncbi`, the module that calls them.
    Patched on the old module, the patch would silently have stopped
    taking effect.
- **Packaging.** `pyproject.toml` lists packages explicitly. Without
  adding `analyses._blast`, `pip install .` would have shipped a
  `blast.py` that cannot import; a clean non-editable install was checked.
- **Rule tests.** The scans in `test_terminal_ui` (layering, UI rules),
  `test_phase6_tier2` (spelling, heuristic vocabulary, output notation)
  and `test_documentation_consistency` read `analyses/*.py` only. They
  now include `analyses/**`, so the moved code stays under the same
  rules. In `test_phase6_tier2` only these scan helpers changed, not the
  CRISPR tests.
- **Ruff.** `analyses/blast.py` ignores F401: its imports are the
  interface.
- **Results.** Unchanged: characterization snapshot, golden master,
  rendered menus and the full suite.

### Package C, step 1: the primer screen without a hidden global (F3)

- **Before.**
  - `_v7_generate_professional_primer_candidates` wrote the design's
    configuration into the module global `_V5_ACTIVE_CONFIG`.
  - `_v5_candidate_quality` read the optimum Tm and GC from it, so a score
    depended on the last design that ran.
  - Two tests set the global from outside.
- **After.** `_v5_candidate_quality(item, config)` takes the
  configuration as an argument. Its two callers, `_v7_screening_score` and
  `_v7_measure_candidates`, pass the `config` they already received; on
  every public path it was the object the global held. The global is
  removed.
- **Evidence.**
  - The characterization snapshot is unchanged: professional design under
    all three Tm models and high salt, restriction-site design, the
    self-tests.
  - The screen-equivalence tests pass.
  - New `tests/test_primer_screen_configuration.py` (3 tests): the
    optimum of the given configuration is used, an earlier design no
    longer changes a later score, and the module keeps no active
    configuration.
- **Tests changed.** `test_primer_corrections` passes the configuration;
  `test_primer_screen_equivalence` no longer sets the global.

### Package C, step 2: `analyses/primer.py` split into `analyses/_primer/`

- **Before.** One module of 5,182 lines in historical layers. The legacy
  engine at the top called functions defined thousands of lines further
  down; Python resolved those calls at run time.
- **After.** `analyses/primer.py` (244 lines) is the public module. It
  keeps its docstring and both `python -m` entry blocks and re-exports
  all 139 names it defines.
- **The modules of `analyses/_primer/`, in dependency order.**
  - `sequence`, `models`, `thermodynamics`, `structure`, `specificity`.
  - `design`: the legacy engine.
  - `professional`: model, scores, single-primer analysis, diagnosis.
  - `screening`: the candidate pool and its Tm-tier screen.
  - `ranking`: pair scoring and diversity selection.
  - `workflow`: `design_primers_professional` and restriction-site design.
  - `selftest`.
  - The import graph has no cycle. Each module reads top-down without a
    forward reference to code elsewhere.
- **How.** With the same script as package B. All 147 statements were
  found verbatim; no code inside a statement changed.
- **The statement that defines no name.**
  - `_SANTALUCIA_2004_STACKS.update({...})` sets the 2004 AA/TT stack. It
    moved into `thermodynamics`, directly after the dict it changes, as in
    the original.
  - A first run of the script had left it in the public module. Through
    the shared dict it would still have taken effect, but only once
    `analyses.primer` was imported.
  - New `InternalModuleIndependenceTests` (in
    `tests/test_module_namespaces.py`) import every internal module of
    both packages without the public module. They check that the 2004
    parameter set is complete. Two further checks: no internal module
    imports a public one, and the public names are the implementation
    objects.
- **Unchanged names.** `GAS_CONSTANT` is still assigned twice with the
  same value. The four `calculate_tm_nearest_neighbor` definitions (three
  overloads and the implementation) stay together.
- **Packaging and lint.** `pyproject.toml` lists `analyses._primer`;
  `analyses/primer.py` ignores F401 like `blast.py`.
- **Results.** Unchanged: characterization snapshot, golden master,
  primer3-py comparisons (run on this machine), screen-equivalence tests
  and the full suite.

### Package G: the CI import check uses the installed package

- **Found.** Every CI step runs in the checkout, so Python imports
  `analyses` from the source tree even after `pip install .`. A package
  missing from the explicit list in `pyproject.toml` would pass every job.
  Package B nearly shipped such a case: `analyses._blast` was not listed.
- **Changed.** The import check in the `tests` job now runs in
  `$RUNNER_TEMP` with `python -I`. It imports `analyses.blast`,
  `analyses.primer` and the earlier modules, asserts that they come from
  `site-packages`, and runs the BLAST self-test there.
- **Checked locally.**
  - The step's commands were run against a clean non-editable install in
    a new virtual environment: they pass.
  - With `analyses/_primer` removed from that installation, the step fails
    with `ModuleNotFoundError: No module named 'analyses._primer'`.
- **Not checked.** The workflow has still not run on GitHub: there is no
  remote.

### Package D: component tests for the least-reached BLAST and primer parts

- **How chosen.** The split showed per responsibility what the suite
  reached least. The tests pin documented behavior; none changes a result.
- **`tests/test_blast_components.py` (19 tests).**
  - Every rejected-input message of query validation and of the FASTA
    loader.
  - Program/query compatibility and the ORF adapter (which CRISPR has its
    own copy of; CRISPR is not touched).
  - The search order for BLAST+ programs and `update_blastdb.pl`.
  - The local engine with a fake `subprocess.run`: parameters to flags,
    `-seg` versus `-dust`, timeout, start failure, non-zero exit, missing
    output, preconditions, and the short-query rescue: it runs once, keeps
    the E-value and states its notes.
- **`tests/test_primer_components.py` (8 tests).**
  - Rejection of invalid reaction conditions and ionic-strength inputs.
  - The specificity screen: arguments, both orientations, both input forms.
  - The professional workflow's warnings with and without references, and
    the failure reason when no pair fits.
- **Observation, not changed.** On a poly-A template the professional
  design returns one pair (A21/T21), flagged `high_risk` and
  `critical_risk` with a heterodimer warning. The screen measures every
  window without hard filters and leaves the risk to the ranking (V6
  section). The test records this as current behavior.
- **Coverage.**

  | Module | Before | After |
  |---|---:|---:|
  | `_blast/query` | 66 % | 97 % |
  | `_blast/executables` | 62 % | 91 % |
  | `_blast/local_engine` | 59 % | 91 % |
  | `_primer/specificity` | 85 % | 100 % |
  | `_primer/thermodynamics` | 85 % | 90 % |
  | `_primer/workflow` | 85 % | 92 % |
  | **Total** | 81 % | 82 % |

- **Suite.** 601 tests, 0 failures on 3.14.6 and 3.10.11.

## Phase 8 – consolidation under the owner's mandate (from 2026-11-18)

The owner's mandate of 2026-11-18 (`docs/PROJECT_HISTORY_AND_DECISIONS.md`
§4) allows Claude, in the copy, to:
- restructure code;
- correct defects that are backed by evidence;
- delete dead code;
- rename internal names;
- drop the ReadOnly protection.

CRISPR and the AI Experiment Planner stay excluded. The rules of
`CONTRIBUTING.md` still hold: the snapshot, the golden master and the
reference tests decide whether a change alters a result.

### Package A: ReadOnly dropped, dead code and legacy files removed

- **ReadOnly.** The flag is cleared on all 42 tracked files that carried
  it. Git history, the test suite and the characterization snapshot are
  the protection; the flag only slowed every edit.
- **Deleted.**
  - `print_protein_report` (`ui/protein_menu.py`, 599 lines) and
    `print_primer_report` (`ui/primer_menu.py`, 39 lines). Neither had a
    caller; they were older versions of reports the menus now draw
    differently.
  - The four imports that only they used.
  - `blast_test.json` (project root): a BLAST run from 2026-10-06 that no
    code or document read.
- **Moved.** `blast_results/TEST_8_REFERENCE_blast.json` and its README
  are now in `tests/fixtures/legacy_blast/` (`SOURCE.md`). The file is
  test evidence: `test_p0_regressions` reads it. The README claimed
  otherwise ("No code or test reads this file"); corrected.
- **`PRIMER_V3_PUBLIC_API`** named two functions that do not exist in
  `analyses.primer`: `print_primer_report` and
  `print_professional_primer_report`, a report printer of
  `ui/primer_menu.py`. Both are removed from the list. A new test checks
  that every listed name exists.
- **Environment, not code.** Windows application control blocked
  ViennaRNA's compiled module on 2026-11-18 ("Eine
  Anwendungssteuerungsrichtlinie hat diese Datei blockiert"). The six
  ViennaRNA tests skip with that reason; the block is not bypassed. They
  passed earlier the same day.

### Package B: the historical `_vN_` primer helpers get descriptive names

- **Before.** The professional primer engine's helpers carried the
  generation they were introduced in (`_v5_`-`_v11_`, `_v3_raw`, `_v4`).
  The public functions were bound to them by assignment
  (`rank_professional_primer_pairs = _v11_rank_professional_primer_pairs`).
  `_V7_NEAREST_NEIGHBOR_METHODS` was a second name for
  `NEAREST_NEIGHBOR_TM_METHODS`.
- **After.** Descriptive names. `generate_professional_primer_candidates`,
  `rank_professional_primer_pairs` and `_primer_ui_design` are defined
  directly; the alias constant is gone.
- **Unchanged.** The "V5"-"V11" section headings, which record history.
- **Mapping.**

| Before | After |
|---|---|
| `_v5_structure_penalty` | `_structure_penalty` |
| `_v5_candidate_quality` | `_screening_quality` |
| `_v6_relaxed_generation_config` | `_measurement_config` |
| `_v7_apply_selected_tm_model` | `_apply_selected_tm_model` |
| `_v7_generate_professional_primer_candidates` | `generate_professional_primer_candidates` |
| `_V7_NEAREST_NEIGHBOR_METHODS` | `NEAREST_NEIGHBOR_TM_METHODS` |
| `_v7_screening_tm` | `_screening_tm` |
| `_v7_screening_score` | `_screening_score` |
| `_v7_screening_upper_bound` | `_screening_upper_bound` |
| `_v7_tm_tier` | `_tm_tier` |
| `_v7_measure_candidates` | `_measure_candidates` |
| `_v7_measure_all_candidates` | `_measure_all_candidates` |
| `_v8_structure_severity` | `_structure_severity` |
| `_v11_model_tm` | `_pair_model_tm` |
| `_v11_tm_optimum` | `_model_tm_optimum` |
| `_v11_piecewise_penalty` | `_piecewise_penalty` |
| `_v11_pair_quality` | `_pair_quality` |
| `_v11_sequence_similarity` | `_sequence_similarity` |
| `_v11_is_diverse` | `_is_diverse` |
| `_v11_select_pairs` | `_select_diverse_pairs` |
| `_v11_rank_professional_primer_pairs` | `rank_professional_primer_pairs` |
| `_v7_primer_ui_design` | `_primer_ui_design` |
| `_primer_ui_design_menu_v4` | `_primer_ui_design_menu` |
| `_calculate_tm_nearest_neighbor_v3_raw` | `_calculate_tm_nearest_neighbor_raw` |

- **Compatibility.** The old names were private. `test_module_namespaces`
  now expects the new names and checks that the old ones are gone.
- **Evidence.** The characterization snapshot (78 scenarios), the golden
  master, the screen-equivalence tests and the suite are unchanged: 603
  tests, 0 failures on 3.14.6 and 3.10.11.

### Package C: `analyses/mutation.py` split into `analyses/_mutation/`

- **Before.** One module of 4,670 lines.
- **After.** `analyses/mutation.py` (307 lines) is the public module.
  It keeps its docstring, `__all__` and the entry block, and re-exports
  all 97 names it defined.
- **The 11 modules of `analyses/_mutation/`.** `sequence`, `models`,
  `alignment`, `detection`, `coding`, `hgvs`, `comparison`, `context`,
  `impact`, `reporting`, `selftest`. No import cycle.
  - `detect_mutations` (alignment, calls, coding consequences and
    statistics in one result) sits in `comparison`, next to
    `compare_sequences`. In `detection` it would have formed a cycle with
    the statistics.
- **How.** Same script and checks as Phase 7.1. All 99 statements were
  found verbatim.
- **Characterization first.** New section `mutation` in
  `tests/characterization.py` with the snapshot
  `tests/fixtures/characterization/mutation.json` (79 scenarios, recorded
  from `727208a` before the split). It covers:
  - five variant types and their reports, JSON, HGVS and alleles;
  - context, restriction and primer impact, regions;
  - variant sets, spectrum, hotspots, statistics, quality;
  - left normalization, alignment, the coding frame, proteins;
  - the advanced analysis and the self-test.
- **Namespace.** `test_module_namespaces` now also guards the 97 names
  and `__all__` of `analyses.mutation`. The independence test includes
  `analyses._mutation`.
- **A test that would have passed vacuously.**
  `test_alignment_memory.test_banded_and_full_alignments_are_identical`
  patched `analyses.mutation.BANDED_ALIGNMENT_MIN_CELLS`. After the move
  that patch no longer reaches the alignment, and the test would have
  compared the banded path with itself. It now patches
  `analyses._mutation.alignment`.
- **Snapshot comparison across interpreters.**
  - CPython 3.10.11 and 3.14.6 on Windows differ in the last digit of
    `math.log2`.
  - `sequence_complexity` of one variant gives 98.82341858018908 on 3.10
    and 98.82341858018907 on 3.14, with the unchanged code of `727208a`
    on both.
  - `tests/test_characterization.py` now compares floats with a relative
    tolerance of 1e-12. Everything else is still compared exactly, and a
    difference of 1e-9 is still reported.
- **Packaging and lint.** `pyproject.toml` lists `analyses._mutation`;
  `analyses/mutation.py` ignores F401.
- **Results.** Unchanged: snapshot, golden master, rendered menus, the
  Mutalyzer HGVS reference test and the full suite. 606 tests, 0 failures
  on 3.14.6 and 3.10.11.

### Package D: `analyses/cloning.py` split into `analyses/_cloning/`

- **Characterization first** (`b9d9d6e`).
  `tests/fixtures/characterization/cloning.json` holds 25 scenarios
  recorded from `8e1f6d2`:
  - planner and best strategies on linear and circular vectors;
  - restriction maps; directional and mismatched strategies;
  - single-enzyme cloning with summary, junctions, reasons and ranking;
  - suitability, compatibility, run modes, overhangs, end compatibility;
  - validation and the self-test.
  `test_module_namespaces` freezes the 61 top-level names. (The commit
  message of `b9d9d6e` says 58; the list has 61.)
- **After.** `analyses/cloning.py` is the public module. The code is in
  seven modules of `analyses/_cloning/`: `cache`, `sites`, `ends`,
  `strategy`, `search`, `summary`, `selftest`. They form a strict layer
  order without a cycle. All 62 statements moved verbatim.
- **The cache stays one module.** `_strategy_search_cache` sets the
  global `_STRATEGY_SEARCH_CACHE` and `_cached` reads it; both are in
  `cache`. `analyses.cloning._STRATEGY_SEARCH_CACHE` is now only the value
  at import time.
- **Three cache tests would have passed vacuously.** `test_planner_cache`
  patched `cloning._analyze_enzyme` and `cloning._cached` and read
  `cloning._STRATEGY_SEARCH_CACHE`. After the move, none of these reaches
  the code: the counting test would have counted nothing and still
  passed. The tests now patch `analyses._cloning.sites` and read
  `analyses._cloning.cache`. The counting test also asserts that it
  counted at least one call.
- **Self-test.** `run_cloning_self_test` checked the public API with
  `globals()`, which were the public module's before the move. It now
  asks `analyses.cloning` (a lazy import inside the function). This is
  the one edit inside moved code; the snapshot's self-test entry is
  unchanged.
- **Re-exported imports.** The golden master and two tests read
  `cloning.RESTRICTION_ENZYMES` and check that `cloning.bases_match` is
  the restriction module's matcher. `analyses/cloning.py` therefore keeps
  importing `RESTRICTION_ENZYMES`, `bases_match` and
  `find_recognition_site_matches`.
- **Packaging and lint.** `pyproject.toml` lists `analyses._cloning`;
  `analyses/cloning.py` ignores F401.
- **Results.** Unchanged: snapshot, golden master, rendered menus and the
  suite. 608 tests, 0 failures on 3.14.6 and 3.10.11; coverage 84 %.

### Package E: distribution name, homepage and the CI import check

- **Distribution name** `synapse` -> `synapse-gen`. The product stays
  **Synapse** (owner decision of 2026-11-18). The name `synapse` is
  already taken on the Python Package Index, and the owner's domain is
  `synapse-gen.com`. Only `pip install` / `pip show` see this name; the
  import names are unchanged. An existing editable install keeps the old
  metadata until it is reinstalled (`pip install -e .`).
- **`[project.urls]`.** Homepage `https://synapse-gen.com`; keywords
  added.
- **CI.** The import check of the installed package now also imports
  `analyses.mutation` and `analyses.cloning` and runs the cloning
  self-test there. It was run locally against a clean non-editable
  install: it passes, and `pip show synapse-gen` reports the homepage.
- **Not done, recorded as the next packaging step.** There is no
  installable `synapse` command yet. `main.py` finds its example files
  next to itself (`examples/`), and they are not part of the installed
  package. A console script needs the examples packaged as data and
  `main` moved into a package.

### Package G: local BLAST left its runtime views behind (defect, fixed)

- **Found** on 2026-11-18, when the owner copied the real `refseq_rna`
  database (21 volumes, 81 GB) into `data/blast_db/`.
- **Defect.**
  - For nested multi-volume databases, `LocalBlastEngine.run` builds a
    runtime alias view of hard links in `synapse_blast_runtime_*` next to
    the database. It never removed it. (`verify_database` does the same in
    a `TemporaryDirectory` and cleans up.)
  - The owner's folder held one from 28 October: 155 hard links and an
    alias file. The database manager then listed it as a second 81-GB
    database, `refseq_rna/synapse_blast_runtime_s3lx6pw3/refseq_rna`.
- **Fix.**
  - `analyses/_blast/local_engine.py` registers the view in an
    `ExitStack`, which removes it when the search ends, also after an
    error. Removing hard links leaves the database files untouched.
  - `analyses/_blast/databases.py` no longer lists anything inside a
    `synapse_blast_*` directory as a database. A view left by an
    interrupted search therefore cannot appear as a database again.
- **Evidence.** New `RuntimeViewCleanupTests` (3 tests in
  `tests/test_blast_components.py`) cover a successful search, a failed
  search and discovery. All three failed on the code before the fix.
- **Owner's data.**
  - The stale view was removed after checking that every file in it was a
    hard link to a file in `refseq_rna.NN/` (the same file ID). The one
    exception was the view's own alias file.
  - `verify_database("refseq_rna")` afterwards: VERIFIED, 83,968,556
    sequences, 249,422,900,906 bases, BLASTDB version 5, dated 30 Sep 2026.
- **Also removed (owner's request).** `data/blast_db/local_test`, a
  byte-identical copy of the tracked `data/local_test`, which the tests
  use. It was sent to the Recycle Bin.

### Package H: an E-value of 0 was read as infinity (scientific defect, fixed)

- **Found** by the first real searches against the owner's `refseq_rna`
  on 2026-11-18. A 400-nt fragment of human GAPDH mRNA (NM_002046.7, nt
  301-700) matched its transcripts with 100 % identity and coverage, yet
  Synapse reported E = infinity and the band `NOT_SIGNIFICANT`.
- **Cause.** `_parse_hsp` read
  `_to_float(Hsp_evalue, math.inf) or math.inf`. BLAST reports E = 0 for
  its strongest matches, and `0 or math.inf` is infinity.
- **Effect, for remote and local searches alike.** Every hit with E = 0:
  - was reported as `NOT_SIGNIFICANT`;
  - was not counted among the significant hits;
  - was ranked behind every weaker hit (`rank_hits` sorts by E-value);
  - was exported as `Infinity`, which is not valid JSON.

  The evidence score itself already treated E <= 0 as the strongest case;
  the parser never passed it on.
- **Fix.** `_evalue_or_infinity` in `analyses/_blast/parsing.py`. A missing
  or unreadable value is still infinite; 0 stays 0.
- **Evidence.**
  - `ZeroEvalueTests` (4 tests): three of them failed before the fix (the
    perfect hit ranked behind a hit with E = 1e-5).
  - The same refseq_rna search after the fix: E = 0.0, band
    `VERY_STRONG_SEQUENCE_EVIDENCE`.
  - Characterization snapshot and golden master unchanged: their XML
    contains no E = 0.

### Package I: opt-in searches against a real refseq_rna

`tests/test_blast_refseq_rna.py` (4 tests) runs Synapse's local engine
against an installed, verified `refseq_rna`, with the GAPDH query above.
It runs only with `SYNAPSE_REFSEQ_RNA_TESTS=1` (about 15 minutes). It
checks:
- the exact fragment: human GAPDH, 100 %/100 %, E = 0, a significant band;
- three substitutions: 99.25 %;
- the reverse complement: a minus-strand hit;
- a 22-nt query: no hit in the standard search, found by the blastn-short
  rescue;
- that no runtime view is left behind.

The CI reference job allows this skip next to the benchmark. Without
that, its gate would fail, which was checked with `check_skips.py`.

Timing on the owner's machine: about 190 s per 400-nt search, and 444 s
for the 22-nt query (standard search plus rescue). blastn used about 7 GB
of memory.

## Phase 9 – completion under the owner's mandate (from 2026-11-19)

The owner's mandate of 2026-11-18 holds: fix all defects and complete the
project; CRISPR and the AI Experiment Planner stay excluded. Windows
Smart App Control was switched off by the owner on 2026-11-19, so
primer3-py, ViennaRNA and `pyright.exe` run again.

**Housekeeping (no code change).**
- The 82-GB BLAST database was moved out of the OneDrive project folder to
  `C:\SynapseData\blast_db`; `SYNAPSE_BLAST_DB_ROOT` (user environment)
  points there. The move was a rename on the same drive; the files were
  local and complete beforehand.
- `refseq_rna` is verified again after the move.
- The stale `refseq_rna/synapse_test.nal` went to the Recycle Bin. It
  pointed to the original "Project Synapse" folder and showed up as an
  extra database.

**Differential checks against Biopython**, on random input
(1,500 sequences for the reverse complements, 400 alignments, 200 GC
values). No defect found:
- all seven reverse-complement implementations equal
  `Bio.Seq.reverse_complement` on IUPAC sequences;
- every mutation alignment reaches the optimal score of
  `Bio.Align.PairwiseAligner` (linear gaps, 2/-1/-2) and preserves both
  sequences;
- GC content is identical;
- the DNA molecular weight differs by the terminal groups (17.7 Da at
  100 nt). That is the documented convention ("terminal-group assumptions
  differ"), not a defect.

### Package A: all 27 NCBI genetic codes

- **Before.** 4 of NCBI's 27 genetic codes (1, 2, 4, 11). The four were
  correct (identical to NCBI and Biopython).
- **After.** All 27, copied from NCBI's gc.prt version 4.6. The file is
  stored verbatim with its SHA-256 in `tests/fixtures/ncbi_gc/`. Code 4
  now carries NCBI's exact name.
- **Context-dependent codons** (codes 27, 28, 31). NCBI lists an amino
  acid for them and marks them as able to terminate; Biopython lists them
  only as stops.
  - `GeneticCode.stop_codons` now includes every codon that can
    terminate, and the new `context_dependent_codons` names these ones.
  - `orf.codon_table_for` ends ORFs at them; `translate` decodes them as
    the amino acid (METHODS).
- **The ORF menu** lists the codes from `available_genetic_codes()`, so
  all 27 are selectable without a UI change.
- **Evidence.** `tests/test_genetic_codes_ncbi.py` (5 tests):
  - the fixture's checksum;
  - every table equal to NCBI's file (amino acids, starts, stops, name);
  - equal to Biopython except exactly the context-dependent codons;
  - ORFs end at TAA under code 28.

  Characterization, golden master and rendered menus are unchanged.

### Package B: ORF finder and RNA translation checked in all 27 codes

- **ORF finder.** `tests/test_orf_differential.py` compares Synapse's ORF
  finder (ATG only, complete ORFs, no nesting) with a deliberately naive
  independent finder: for each stop, the most upstream in-frame ATG after
  the previous stop, on both strands.
  - On 2026-11-19 the two agreed on all 675 random sequences (25 per code,
    60-900 nt, minimum lengths 30, 75 and 150).
  - The suite runs 4 per code plus a fixed case. That case guards against
    agreeing on empty sets.
- **RNA translation.** `rna.translate_rna_frame` equals the NCBI tables up
  to the first stop in 1,620 random frames, across all codes and all three
  frames.
- **Biopython agrees on the context-dependent codons.** `Bio.Seq.translate`
  decodes them as amino acids too, as `GeneticCode.translate` does: "Such
  codons will be translated as amino acid".
- No defect found.

### Package C: size, release and disk space before an NCBI database download

- **Before.** The database manager started `update_blastdb.pl` after the
  user typed DOWNLOAD. It showed no size and did not check the disk:
  core_nt needs about 308 GB, refseq_protein about 253 GB.
- **After.** New `analyses/_blast/ncbi_metadata.py` reads NCBI's metadata
  file for the chosen database (a few kB; JSON, metadata version 1.1).
  - Before the confirmation, the dialog shows NCBI's release date, the
    number of sequences, the download size, the size on disk, the free
    space needed and the free space on the drive.
  - It says whether NCBI's release is newer than an installed copy (the
    date `blastdbcmd -info` reports).
  - If the drive is too small, no download is offered, and the message
    points to `SYNAPSE_BLAST_DB_ROOT`.
  - `update_ncbi_database(..., required_free_bytes=...)` refuses to start
    as well.
  - Without network access the dialog says the size is unknown and still
    allows the download, as before.
- **Space estimate.** The decompressed size, less what an installed copy
  occupies, plus one average volume archive: `update_blastdb.pl` unpacks
  one volume at a time.
- **Owner's data on 2026-11-19.** NCBI's refseq_rna release is 2026-10-30
  (85,490,457 sequences); the installed copy is from 2026-09-30. The
  dialog would offer the update.
- **Evidence.** `tests/test_blast_download_plan.py` (8 tests) with NCBI's
  real swissprot metadata file as a fixture (SOURCE.md).
- **Interface.** The new names are exported by `analyses.blast`; the
  menu imports them from there.

### Package D: an installable `synapse` command, and no user files in site-packages

- **Before.** `pip install .` installed the packages but no way to start
  the session: `main.py` was not packaged, and the example FASTA lived
  outside any package. An installed Synapse would also have written into
  `site-packages`:
  - logs to `site-packages/logs`;
  - session exports to `site-packages/exports`;
  - the default BLAST database root, `site-packages/data/blast_db`
    (`_resolve_default_local_db_root` found the installed `data`
    package).
- **New `core/paths.py`.**
  - From the source folder (detected by `pyproject.toml` next to the
    packages), Synapse's folder is the project folder, as before.
  - Installed, it is a per-user folder: `%LOCALAPPDATA%\Synapse`, or
    `$XDG_DATA_HOME/synapse` (`~/.local/share/synapse`) elsewhere.
  - `SYNAPSE_HOME` overrides it; `SYNAPSE_LOG_DIR` and
    `SYNAPSE_BLAST_DB_ROOT` still win for their own files.
  - Logs, exports, relative paths and the BLAST root use it.
- **The session moved.** `main.py` -> `synapse_gen/app.py` (`git mv`).
  - `main.py` is now a launcher: `import main` returns the app module
    itself, so `main.main()` and patches of `main` keep working.
  - New `cli()` wraps `main()`. When the input ends (Ctrl+Z / Ctrl+D or
    the end of piped input) or on Ctrl+C, the session now ends with a note
    instead of a traceback (exit code 0, or 130 for Ctrl+C).
  - `[project.scripts] synapse = "synapse_gen.app:cli"`, plus
    `python -m synapse_gen`.
- **The example moved.** `examples/dna.fasta` -> `data/examples/dna.fasta`,
  installed as package data. `examples/crispr_batch_test.fasta` stays,
  because CRISPR is excluded.
- **Found while doing this.**
  1. `_test_default_database_root_resolution` only accepted a root ending
     in `data/blast_db`. With the owner's `SYNAPSE_BLAST_DB_ROOT` (set on
     2026-11-18), the BLAST self-test, and with it the suite, would have
     failed in every new terminal. It now checks the exact expected root
     for each configuration, and that the root is never inside the
     `analyses` package.
  2. The suite depended on the user's environment. `tests/__init__.py` now
     removes `SYNAPSE_BLAST_DB_ROOT`, `SYNAPSE_HOME`, `SYNAPSE_LOG_DIR` and
     the NCBI identification variables. It keeps them in
     `USER_ENVIRONMENT` for the opt-in refseq_rna test, which needs the
     real database folder.
  3. My own Phase 8 change to the CI import check put
     `["all_passed"]` inside a double-quoted `python -c "..."`, which ends
     the shell string. It was corrected to `['all_passed']`. The step is
     now verified by running the workflow's lines exactly as written
     (extracted from the YAML) against a clean installation.
- **Tests.**
  - `tests/test_app_entry.py` (4): the launcher, the clean ends, the
    declared command.
  - `StorageLocationTests` (5, in `test_blast_components`): source folder,
    installed, precedence of the settings, logs, the packaged example.
  - The tests that scan or start `main.py` now read
    `synapse_gen/app.py` where they check the session's code.
- **Verified on a fresh installation.** In a new virtual environment with
  a non-editable install:
  - `synapse` starts and loads the bundled example;
  - logs land in `SYNAPSE_HOME/logs`;
  - the CI import step passes.
- **Not changed.** CRISPR exports still default to the project folder
  (`analyses/crispr.py`, excluded); set `SYNAPSE_EXPORT_DIR` when
  Synapse is installed.

### Package E: `analyses/protein.py` split into `analyses/_protein/`

- **Characterization first.** The section `protein`
  (`tests/fixtures/characterization/protein.json`) holds 103 scenarios
  recorded from `ee8203d`:
  - four proteins (one with ambiguous residues and a stop), with every
    property function, the combined analyses and the charge at four pH
    values;
  - four RNAs, among them IUPAC codons, translated under codes 1, 2 and 11;
  - codon resolution.

  `test_module_namespaces` freezes the 56 names.
- **After.** Seven modules: `constants`, `sequence`, `translation`,
  `composition`, `charge`, `profiles`, `report`. No cycle; all 56
  statements moved verbatim. `analyze_protein`, the composition count,
  sits in `composition`; in `report` it formed a cycle.
- **Snapshot comparison.** Floats now also allow an absolute difference
  of 1e-12. Since Python 3.12, `sum()` of floats is compensated, so a
  hydropathy window that sums to 0 gives -7.4e-17 on 3.10 and -8.6e-17
  on 3.14.
- **Plan.** `scripts/refactor/plans/protein.py`. Packaged in
  `pyproject.toml`; covered by the independence test.

### Package F: `analyses/variant_analysis.py` split into `analyses/_variant/`

- **Characterization first.** The section `variant`
  (`tests/fixtures/characterization/variant.json`) holds 45 scenarios
  recorded from `6e52ca7`. It runs four VCFs (the two module samples and
  two of the bcftools test vectors) through:
  - analysis, report, VCF and TSV export;
  - statistics, cohort summary, prioritization;
  - normalization with contig references in both matching modes;
  - a quality filter, set comparison, intervals, feature overlap;
  - the self-test.

  The 84 names are frozen.
- **After.** Eight modules: `models`, `vcf`, `normalization`,
  `annotation`, `statistics`, `export`, `workflow`, `selftest`. The
  public module keeps `__all__`. No cycle; all 86 statements moved
  verbatim.
  - `classify_variant_type` sits in `models`, because `Variant` calls it.
  - `VariantRecord` sits in `vcf`, beside the remapping helpers it uses.
- **Results.** Unchanged on 3.14.6 and 3.10.11, including the bcftools
  vectors.

### Package G: sequence complexity depended on Python's hash seed (defect, fixed)

- **Found** while recording the characterization snapshot of
  `analyses/sequence.py`: two runs of the same code gave different
  results.
- **Cause.** `calculate_sequence_complexity` summed the Shannon entropy
  terms while iterating over the set `{"A", "T", "C", "G"}`. Python
  randomizes set iteration per process (`PYTHONHASHSEED`), so the order
  of the four additions changed between runs, and with it the last digit
  of the complexity.
- **Visible effect.** Where windows tie, the reported position of the
  complexity maximum depends on that last digit; one test sequence gave
  position 180 in one run and 145 in the next. The overview report was
  therefore not reproducible.
- **Fix.** The bases are iterated in sorted order. The values change at
  most in the last digit.
- **Evidence.**
  - New `tests/test_determinism.py` runs a sequence analysis under four
    hash seeds and requires identical output. It failed before the fix.
  - The same comparison of all 7 characterization sections under three
    hash seeds found no other dependence; only the sequence section
    differed (252 values).
  - The characterization section `sequence` (111 scenarios: five sequences
    including IUPAC, CpG-rich and homopolymer cases, every overview
    function, the profiles, the self-test) was recorded after the fix. It
    is identical under four further hash seeds and on 3.10.11.

### Package H: `analyses/sequence.py` split into `analyses/_sequence/`

- **After.** Six modules in strict layer order: `bases`, `composition`,
  `profiles`, `cpg`, `summary`, `selftest`. All 43 statements moved
  verbatim; the 42 names are frozen.
- **Results.** The snapshot recorded after package G is unchanged; so are
  the golden master and the rendered menus. Suite OK on 3.14.6 and
  3.10.11.

### Package I: `analyses/orf.py` split into `analyses/_orf/`

- **Characterization first.** The section `orf`
  (`tests/fixtures/characterization/orf.json`) holds 27 scenarios
  recorded from `b46ed3d`:
  - pUC19, linear and circular, under five settings: default, legacy,
    code 11 with table starts, nested ORFs under code 2, 5'-partial ORFs;
  - the enhanced analysis, relationships, ranking, comparison, the feature
    summary, selected ORFs and their translations;
  - an IUPAC case, codon tables of codes 1, 4 and 28, start-codon sets.

  It is identical under two hash seeds and on 3.10.11. The 54 names are
  frozen.
- **After.** Five modules: `codes`, `features`, `finder`, `summary`,
  `selection`. No cycle; all 54 statements moved verbatim. One generated
  import was removed: `finder` imports `DNA_IUPAC` from `core.alphabets`
  inside a function.
- **`scripts/refactor/split_module.py`** no longer counts names a
  function imports itself as module-level uses.
- **Independence test.** It now checks that an internal module never
  imports its own public module. Other analyses are allowed: `_orf` uses
  `analyses.sequence`.
- **Results.** Unchanged: ORFfinder comparison, the differential ORF test
  in all 27 codes, the golden master.

### Package J: `analyses/restriction.py` split into `analyses/_restriction/`

- **Characterization first.** The section `restriction`
  (`tests/fixtures/characterization/restriction.json`) holds 46 scenarios
  recorded from `709f964`. On pUC19, linear and circular, with 10 enzymes,
  it covers:
  - the analysis, maps, summaries, digests, gel bands, cut boundaries;
  - unique and non-cutters, pair recommendations, double digests;
  - every enzyme's sites;
  - an IUPAC digest with possible sites, the enzyme-database check,
    cloning-compatible sites and the self-test.

  It is identical under two hash seeds and on 3.10.11. The 51 names are
  frozen.
- **After.** Five modules: `sites`, `digest`, `compatibility`, `api`,
  `selftest`. No cycle; all 52 statements moved verbatim.
- **Two adjustments, the same as for cloning.**
  - `run_restriction_self_test` checked its public API with `globals()`;
    it now asks `analyses.restriction`.
  - `RESTRICTION_ENZYMES`, `get_restriction_enzyme` and
    `get_type_iis_enzyme` stay importable from `analyses.restriction`.
    The restriction menu imports them there; 31 tests failed until they
    were re-exported.

### Package K: `analyses/rna.py` split into `analyses/_rna/`

- **Characterization first.** The section `rna`
  (`tests/fixtures/characterization/rna.json`) holds 101 scenarios
  recorded from `8dff342`: four sequences (one with IUPAC symbols, one
  hairpin-rich, one with homopolymers) under codes 1 and 2, every
  analysis function, the profiles and the pipeline source.
- **Determinism.** The exact comparison under three hash seeds now covers
  all ten sections and finds no difference.
- **After.** Six modules: `bases`, `composition`, `coding`, `structure`,
  `summary`, `pipeline`. `_orf_value` moved to `bases`, beside its caller
  `orf_genetic_code`, to avoid a cycle. All 47 statements moved verbatim;
  the 47 names are frozen.
- **Results.** Suite OK on 3.14.6 and 3.10.11; the ViennaRNA comparisons
  ran.

### Package L: the CI imports every installed module; `literature.py` stays one module

- **CI.** The import check of the installed package walks every module of
  `analyses`, `core`, `data`, `fileio`, `synapse_gen` and `ui` (138 on
  2026-11-19) instead of a hand-kept list, so a package missing from
  `pyproject.toml` fails the job. It also runs the restriction self-test.
  The step was run exactly as written, against a clean non-editable
  install.
- **`analyses/literature.py` is not split, on purpose.**
  - Its tests replace `_request_json`, `_request_text`, the provider
    searches, `urllib.request.urlopen` and `time.sleep` on the module 20
    times.
  - The providers share one module-level throttle (`_LAST_REQUEST_AT`,
    set with `global`).
  - After a split, a replacement on the public module would no longer
    reach the code that calls it, while the tests would still pass. That
    is the failure this phase found three times.
  - The module is cohesive: one search pipeline over four providers. The
    expected gain does not justify the risk.

## Phase 10 – genome scale, curated benchmarks, release preparation (2026-11-20)

The owner's mandate of 2026-11-18 holds (fix every defect, complete the
project; CRISPR and the AI Experiment Planner stay excluded; the licence
and the GitHub repository are the owner's). Seven commits from `6d41d8a`
to `5b1dc8b`, then the documentation. Measurements are in
`docs/BENCHMARKS.md`.

### A. Genome-length input (`6d41d8a`)

*E. coli* K-12 MG1655 (NC_000913.3, 4,641,652 bp) through every analysis
(`scripts/benchmarks/genome_scale.py`). Four analyses did not scale:

- **ORF search, 499 s.** `analyze_orf_relationships` compared all
  n(n−1)/2 ORF pairs (32,406 ORFs). An interval sweep visits only ORFs
  whose start can overlap, in the original order.
- **Restriction digest with all 22 enzymes, more than 15 minutes.** Each
  fragment scanned every cut for its two boundaries; the cuts are now
  indexed by coordinate. Sites are found with a regular expression built
  from exactly the symbols `bases_match` accepts; the sequence check
  enumerates positions only when an unusual symbol occurs.
- **Sequence overview, 154 s and 3.5 GB.** CpG islands and the GC and
  complexity profiles recounted each window; they use prefix counts
  (`core/windows.py`). Above 200,000 windows the stored profiles keep
  every k-th window (`profile_resolution` says so); the statistics still
  use every window.
- **RNA analysis, more than 10 minutes.** The hairpin screen extends each
  stem outward from its innermost pair once and stops at the first
  mismatch no candidate can exceed; the G-U figure is counted without a
  list; the RNA profiles use prefix counts and the same thinning.

Every change returns the former result: differential tests against the
former loops (`test_window_scans`, `test_rna_window_scans`,
`test_orf_relationship_sweep`, `test_restriction_fragment_index`), and an
end-to-end comparison with the code of `ca6d941` on *E. coli* segments
(overview, three ORF settings, RNA, two digests): bit-identical results,
3-23 times faster. Behavior change for long input only: the thinned
stored profiles (documented in `docs/METHODS.md`).

### B. HGVS benchmark on ClinVar; allele fix (`1fd914f`)

1,500 ClinVar variant names on six RefSeq transcripts
(`tests/fixtures/clinvar/`, public domain, SHA-256 checked). c.: 1,489
identical, 11 split delins equivalent; p.: 1,186/1,186.

- **Defect.** `hgvs_c_allele` shifted every change of a multi-change
  allele to its most 3′ position as if it were alone. For PTEN
  c.528_534delinsA it wrote `c.[532_534del;531_532del;534del]`,
  overlapping deletions that describe no sequence. The shift now stops
  one unchanged base before the next change (`shift_limit`). Single
  changes are unchanged.

### C. Per-ORF features (`f6a5d92`)

Codon usage is counted once per candidate without a per-codon loop, the
diversity reuses it, the ambiguity check returns at once for plain ACGT,
and the scan caches codon status. Genome, default settings: 35 s → 11 s;
alternative starts at 90 nt: 276 s → 87 s. Same ORFs.

### D. Dialog tests (`8696b1f`)

16 scenarios for the protein, restriction and cloning dialogs. Two
defects: three protein errors lacked the ERROR prefix; Back inside a
cloning analysis left the cloning menu.

### E. Literature authors (`713a1d4`)

Europe PMC's "Smith A" did not match PubMed's "Ann Smith", so merged
records listed authors twice in results and exports. A trailing token of
one to three capital initials is now also read as initials.

### F. `analyses/literature.py` split (`4bb71d0`, `5b1dc8b`)

Snapshot with the network replaced at `urlopen`, then ten modules in
`analyses/_literature/`; stand-ins moved to the calling modules. The
splitting tool kept only one of several `import urllib.x` lines; fixed.

### Benchmarks on curated data

- ORFs against the RefSeq annotation of *E. coli* K-12: 4,250 of 4,294
  CDS found at their stop codon; every CDS of at least 90 nt
  (`scripts/benchmarks/orf_refseq.py`).
- Primer Tm against published measured values: not done. No measured
  dataset with a clear licence was obtained, and none was reconstructed
  from memory. The comparison with Primer3 stands.

### Documentation and release preparation

New: `docs/USER_GUIDE.md` (checked by `tests/test_user_guide_examples.py`),
`docs/BENCHMARKS.md`, `docs/RELEASE_NOTES.md`, `CITATION.cff` (valid CFF
1.2.0), issue and pull-request templates, Dependabot for the workflow
actions, `data/examples/pUC19.fasta`. Updated: README, METHODS,
FEATURE_STATUS, ARCHITECTURE, TOOL_AND_DATA_LEDGER, TEST_COVERAGE,
PROJECT_HISTORY_AND_DECISIONS, `scripts/refactor/README.md`.

Correction to the commit message of `6d41d8a`: it says the RNA analysis
became "about 7x faster on 300 kb". Measured without a profiler, a
300-kb RNA analysis took 29.4 s after the early stop and 10.7 s after the
outward extension; the time of the unchanged code was not measured
without a profiler, so the factor is not established. The genome figure
in `docs/BENCHMARKS.md` (more than 600 s before, 137 s after) is measured.
