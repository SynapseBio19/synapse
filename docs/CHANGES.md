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

- Git repository initialised; tag `baseline-audit` = project exactly as audited.
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
  a 99-bp query. Evidence bands ignored significance: E = 7.45 was labelled
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
  are not recognised as starts by the mutation module.

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

### D5 Cloning analysed circular vectors as linear (`analyses/cloning.py`, `core/pipeline.py`, `analyses/restriction.py`)

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
  this commit changes no existing behaviour.
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
- Behaviour change, intended: DNA records with IUPAC ambiguity codes
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

### D7 The stale BLAST artifact is labelled as legacy (`blast_results/README.md`)

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

### D8 CRISPR guides report the modelled cut site (`analyses/crispr.py`)

- Problem: no guide carried a cut position; the only cleavage information
  was a free-text `cut_model` string.
- Now: each system defines where it cuts, counted along the guide in
  targeting orientation, with its source: SpCas9 and SpCas9-NG blunt
  after guide position 17 (3 nt 5' of the PAM; Jinek et al. 2012), SaCas9
  blunt after position 18 (3 nt 5' of the PAM, modelled as for SpCas9;
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
  behaviour (window 4–8 = the CBE/ABE definitions; GC ranked; low-quality
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

### D15 The unimplemented planner is labelled as such (`main.py`, `ui/menu.py`)

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
