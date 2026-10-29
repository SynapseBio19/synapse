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
