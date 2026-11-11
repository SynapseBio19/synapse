# Synapse – methods, defaults and limitations

This document states what each analysis computes, which defaults it uses,
where those defaults come from, and what has (and has not) been checked.
Synapse is research and teaching software. **No result has been validated
experimentally, and Synapse has not been benchmarked against established
tools except where a specific check is listed below.** Context-dependent
thresholds are conventions, not universal biological limits; all of them can
be changed.

Validation levels used below:

- **Known-answer test**: reproduces a published worked example or a value
  computed by hand from the published equations.
- **Implementation cross-check**: agrees with an independent implementation
  of the same equations written for the test suite (checks the code, not the
  model).
- **Not independently validated**: covered only by regression (golden) tests,
  which detect changes but do not prove correctness.

Coordinates: user-facing positions are 1-based and inclusive. Internally,
`core.models.Interval` is 0-based half-open with a strand; Synapse ORF dicts
use 0-based inclusive `start <= stop` on the forward strand.

---

## Sequence input (`fileio/`)

- FASTA and GenBank (`.gb/.gbk`). GenBank topology is read from the LOCUS
  line; FASTA headers may carry `topology=circular`.
- Molecule type: GenBank LOCUS (`aa` = protein; `DNA`, `mRNA`, `ss-RNA`, …
  = nucleotide) and the NCBI FASTA modifier `[moltype=...]`
  (`[moltype=protein]` is a Synapse extension) are read as a declared
  type. Each record is classified with `core.alphabets` and carries its
  type and its evidence level into provenance and the report. Without a declared
  type the letters decide (see "Sequence types").
- Main menu: a record whose type cannot be determined is not guessed; the
  user declares DNA, RNA or protein. RNA records are analyzed in the DNA
  alphabet (U read as T), which the report records. Protein records open
  the protein menu (protein analysis, BLAST, report).
- `core.api` entry points that make external requests refuse to run
  unless the caller confirms explicitly: `blast_search(engine="remote")`
  needs `send_to_ncbi=True` and `literature_search` needs
  `send_to_providers=True`. This keeps the privacy gate the terminal menu
  applies, which a non-interactive layer cannot ask for.
- `core.api`: nucleotide analyses accept DNA and RNA (RNA converted, with a
  warning); restriction digests and CRISPR guide design accept DNA only,
  because both act on DNA; `analyze_protein` takes proteins. Every entry
  point accepts `sequence_type=` for an undetermined sequence.
- Provenance: file SHA-256, sequence SHA-256, path, format, record index.
- The GenBank reader reports what it could not read as record notes:
  symbols dropped from ORIGIN (S21) and feature locations in forms it
  does not parse, such as `order()`, remote accession references and
  between-base sites, which keep their raw string and have no interval
  (S10).
- GenBank locations: simple spans, `complement()` and `join()` are parsed;
  other forms keep only the raw string.

## Sequence types (`core/alphabets.py`)

- Alphabets: IUPAC nucleotide codes (Cornish-Bowden 1985) with T for DNA
  and U for RNA; the IUPAC-IUB one-letter amino-acid codes (JCBN 1984)
  with the ambiguity codes B, Z, X and the later additions J (Leu/Ile),
  U (selenocysteine) and O (pyrrolysine), as used by NCBI and Biopython;
  `*` for a stop. Gaps and digits belong to no alphabet. The letter sets
  are checked against `Bio.Data.IUPACData` in the tests.
- Normalization removes whitespace and uppercases; it never converts U to
  T and never drops a symbol.
- Classification without a declared type (`classify_sequence`) gives a
  sequence type and an evidence level (`evidence`: CONFIRMED, PROBABLE,
  AMBIGUOUS or INVALID). The level is categorical; it is not a
  probability:
  - symbols that exist only in the protein alphabet (E F I J L O P Q X Z
    `*`) → protein, CONFIRMED;
  - only A, C, G, T/U and N → DNA or RNA, PROBABLE (never confirmed:
    every nucleotide code is also an amino-acid letter);
  - other nucleotide ambiguity codes, or T and U together → AMBIGUOUS:
    no type is chosen and the caller must declare one;
  - empty, or symbols outside every alphabet → INVALID.
- This is a categorical rule about which letters occur, not a cut-off on
  their frequency. No published threshold exists for telling a nucleotide
  sequence from a peptide by letter frequency, so Synapse uses none.
- A declared type (by the caller, or by file metadata) is CONFIRMED when
  the sequence fits that alphabet and INVALID when it does not. File
  metadata usually says only "nucleotide": GenBank and other databases
  write RNA molecules in the DNA alphabet, so the letters then decide
  between DNA and RNA.

## Sequence overview (`analyses/sequence.py`)

- GC and AT content are percentages of the unambiguous bases. With IUPAC
  ambiguity codes present, a minimum and maximum over every resolution of
  the codes are reported as well.
- GC skew = (G − C)/(G + C); AT skew = (A − T)/(A + T); 0 when the
  denominator is 0.
- Sequence complexity is the Shannon entropy of the A/C/G/T composition
  divided by its maximum of 2 bits, so it runs from 0 (one base) to 1
  (all four equally frequent). It measures composition only and does not
  detect repeats: `ACGT` repeated 100 times scores 1.000 and `AT`
  repeated 100 times 0.500. The classes high (≥ 0.75), moderate (≥ 0.50)
  and low are Synapse conventions. The RNA, CRISPR and mutation modules
  use the same measure (the mutation module as a percentage); the field
  names keep the word "complexity" for compatibility.
- Homopolymers: runs of at least 4 identical bases.
- CpG islands: 200-bp windows with GC ≥ 50% and an observed/expected CpG
  ratio ≥ 0.60, where expected = C·G/window length (the criteria of
  Gardiner-Garden and Frommer 1987).
- Molecular mass: composition-based sums of nucleotide residue masses
  for single- and double-stranded DNA, with a range for ambiguous bases;
  the average estimate uses 660 Da per bp. Calculators differ in their
  residue-mass and terminal-group conventions, so the value is an
  estimate.
- Validation: the module self-test (`run_sequence_self_test`, run by the
  suite) checks counts, content, skew, entropy, CpG count, IUPAC ranges
  and the average mass on small constructed sequences, and the golden
  master guards the outputs. CpG-island windows and the residue masses
  are not compared with an external tool (not independently validated).

## ORF search (`analyses/orf.py`)

| Parameter | Default | Alternatives |
|---|---|---|
| Minimum length | 75 nt incl. stop | any value ≥ 3 |
| Start codons | ATG only | ATG/GTG/TTG (`atg_and_alternative`); NCBI starts of the table (`table_starts`) |
| Nested ORFs | longest ORF per stop codon | every in-frame start |
| Ranking | length | legacy heuristic score |
| Genetic code | NCBI table 1 | tables 2, 4, 11 |

- Defaults follow the conventions of NCBI ORFfinder (minimal length 75 nt,
  ATG only, outermost start per stop). They are screening conventions:
  short ORFs, non-ATG starts and nested ORFs exist and can be enabled.
- `LEGACY_ORF_SETTINGS` reproduces the pre-2026-10-29 behavior (10 nt,
  ATG/GTG/TTG, all starts, heuristic ranking); a test proves equivalence.
- The first codon of an ORF is translated as initiator Met when it is an
  initiation codon of the ORF's genetic code (`orf.initiation_codons`:
  the NCBI start codons of the table plus, for table 1, Synapse's GTG
  policy; table 1 = ATG, CTG, GTG, TTG). An IUPAC codon counts only if
  every reading is an initiation codon. Internal codons always use the
  elongation code (e.g. internal CTG = Leu, internal GTG = Val).
- An ORF is translated with its own genetic code: `translate_orf_sequence`
  and `core.api.analyze_protein_from_cds` accept the ORF dict and read its
  `genetic_code` (a conflicting explicit table is an error). A bare DNA
  string defaults to table 1.
- The "heuristic score" (`orf["heuristic_score"]`, 0-100 with a
  High/Moderate/Low/Very Low band) is a ranking aid, not a probability of
  coding. It adds fixed points for completeness, start codon, length,
  translation safety, ambiguous codons, Kozak context and a
  Shine-Dalgarno match; the weights are Synapse conventions.
- Partial ORFs (P1-9, 2026-11-11). 3′-partial ORFs (a start codon but no
  stop before the end) are reported by default. 5′-partial ORFs, open frames
  entering from the 5′ end whose start codon lies upstream of the sequence,
  are reported on request ("Partial - Missing Start", `include_5prime_partial`;
  in the menu the Partial ORFs setting cycles 3′ end / 5′ and 3′ ends /
  excluded). Their first codon is not a start codon, so it is never read as
  initiator Met (ORF translation, protein analysis, `core.api`). Complete ORFs
  inside such a frame are kept; the per-stop rule compares ORFs of one
  status.
- Circular sequences (P1-9, 2026-11-11). The search is repeated on
  dna + dna[:n−1], and an ORF that starts in the first copy and ends in the
  second is reported as running across the origin (`wraps_origin`). Its
  coordinates stay unrolled (start ≤ stop, so length and sequence hold);
  `stop_on_circle` and the displays give its last base on the circle. Partial
  ORFs are not reported on a circle, where they would be artifacts of the
  arbitrary origin. The pipeline and `core.api` pass the record's topology.
  The longest-ORF-per-stop rule keys on strand, 3′ end and status; the frame
  was dropped from the key because it changes across the origin of a circle
  whose length is not a multiple of 3 (on a linear sequence strand and 3′
  end fix the frame, so linear results are unchanged).
- Validation: known-answer tests for frames, strands, coordinates and the
  table-2 read-through of TGA; for circles, an ORF across the origin on both
  strands and rotation invariance (the same ORFs for four origin positions of
  a 900 bp circle); not compared against ORFfinder output (no ORFfinder
  binary; its web service would need an external request).
- Organism mode (P1-9) is not implemented as a separate switch: the genetic
  code and the start-codon policy cover the cases (e.g. table 11 with
  `table_starts` for bacteria), and a preset would only bundle them.

## RNA and protein (`analyses/rna.py`, `analyses/protein.py`)

- Translation in frame +1; ambiguous codons give `?` unless all expansions
  encode one amino acid. With `initiator_start=True` (protein analysis of
  a selected ORF, `core.api.analyze_protein_from_cds`) the first codon is
  decoded as Met only if it is an initiation codon of the genetic code —
  the same rule as the ORF translator.
- Protein analysis uses the genetic code of the selected ORF, and so
  does the RNA module's codon view (see below).
- Physicochemical values (mass, pI, GRAVY and similar) are sequence-based
  estimates. Checked against ExPASy ProtParam for ubiquitin (MW 8564.84 Da,
  pI 6.56, GRAVY −0.489, ε 1490 M⁻¹cm⁻¹).
- External reference comparison (2026-11-02, Biopython 1.88
  `Bio.SeqUtils.ProtParam`, 31 proteins of 41-400 aa): pI within 1e-4,
  GRAVY and aromaticity identical, extinction coefficients identical,
  instability index within 2e-4, average mass within 2e-5 relative (the
  residue mass tables differ slightly). `tests/test_external_reference.py`
  repeats it, skipped without Biopython.
- Instability index (Guruprasad et al. 1990): 10/L times the sum of the
  dipeptide weights. The weight for the dipeptide IG was −7.49 in Synapse
  and 1.0 in Biopython's table (`ProtParamData.DIWV`, which cites
  Guruprasad et al. 1990); it is now 1.0 (2026-11-02). The only remaining
  difference is FY, 33.6 here and 33.601 there. The primary table was not consulted for this correction
  (no external access in this session), so the reference is one
  secondary source, not the paper; confirming it is open work.
- Unresolved residues ("?", from ambiguous codons) keep their position.
  Sliding windows (hydropathy, flexibility, transmembrane screen, signal
  peptide) and the instability index skip the windows or dipeptides that
  contain one instead of deleting the residue, so reported positions always
  refer to the real protein. Composition-based values (pI, GRAVY, extinction
  coefficient and similar) are computed over the resolved residues only.
- The transmembrane and signal-peptide outputs are screens based on the
  Kyte–Doolittle scale, not predictions from TMHMM/Phobius or SignalP.
- RNA hairpin screen: stems of 4-7 nt around loops of 3-12 nt are reported
  when at least 75 % of the stem positions form Watson-Crick pairs. It is a
  sequence screen without free energies (no Turner parameters, no RNAfold).
  G·U wobble pairs (P2-5, 2026-11-11) are counted in a second figure,
  "Including G-U pairs" (`hairpin_count_with_gu`, `allow_gu=True`), not in
  the screen itself: treating G·U like G·C without an energy model raises
  the candidates on the golden-master reference from 58 to 260 rather than
  finding stronger stems.
- RNA codon view (P1-9): the three-frame preview uses the genetic code of
  the selected ORF, taken from the ORF itself; without a selected ORF it
  uses table 1. The +1 frame therefore agrees with the ORF's own
  translation.
- Protein letters (D14-F6): every position is kept. B, Z, J, X and `?`
  are ambiguous; U (selenocysteine) and O (pyrrolysine) are defined
  residues with masses (Biopython `IUPACData.protein_weights` minus
  water); one terminal `*` is a stop. Position-dependent scales skip
  windows containing any non-standard residue. Molecular weight is a range
  for B/Z/J and undetermined for X and `?`. Hydropathy, charge and pI have
  no values for U and O.

## Restriction and cloning (`analyses/restriction.py`, `analyses/cloning.py`)

- Sites are CONFIRMED (exact bases) or POSSIBLE (match only through IUPAC
  ambiguity). Digests use CONFIRMED sites by default
  (`include_possible=False`).
- Restriction and cloning use one matcher
  (`restriction.find_recognition_site_matches`). Only CONFIRMED sites define
  cut positions and cloning ends. An enzyme with a POSSIBLE site in a
  molecule is treated as neither a confirmed unique cutter nor a confirmed
  non-cutter there (the number of real sites is undetermined), so it cannot
  define a cloning end in that molecule.
- Cloning takes the topology of vector and insert (`vector_topology`,
  `insert_topology`; the pipeline and menus use the loaded record's).
  On a circular molecule the segment between the two selected sites, in
  coordinate order, is taken as the region replaced (vector) or excised
  (insert); the result states this assumption.
- Circular topology finds sites spanning the origin, in the digest and in
  every site-based result (map, per-enzyme results, unique and non-cutters,
  enzyme pairs). Positions on a circular sequence are reported modulo its
  length; a cut exactly at the origin is position 0 (1-based 1).
- Coordinates in site and cut records: `site_position`, `top_cut` and
  `bottom_cut` are 1-based (for a cut, the first base after it); the
  `*_0based` keys hold the 0-based site start and the cut boundary (the
  number of bases before the cut). `cut_positions` and fragment boundaries
  are cut boundaries. Until 2026-11-11 digest events used the unsuffixed
  names for 0-based values while maps used them for 1-based values; both now
  follow the map (golden master unchanged, as it stores map records and cut
  boundaries).
- Type IIS enzymes (BsaI, BsmBI, BbsI, SapI; `TYPE_IIS_ENZYMES`) cut outside
  their asymmetric site. The digest searches the site and its reverse
  complement; for a forward site at p the cuts are p + cut_top and
  p + cut_bottom, for a reverse site they are mirrored to the other side.
  The overhang is the top-strand sequence between the cuts, so it depends on
  the molecule; for that reason Type IIS enzymes are not part of the cloning
  compatibility model, which assigns each enzyme one fixed overhang, and
  Golden Gate assembly planning is not implemented. On a linear molecule a
  site whose cuts would fall beyond an end does not cut.
- End compatibility requires the same overhang sequence *and* polarity
  (5′ vs 3′); blunt ends are compatible with blunt ends.
- Methylation sensitivity is shown as a note (ClaI, XbaI, DpnI, HpaII);
  blocking by methylation and star activity are not modeled.
- Validation: known-answer tests (ClaI cut positions and overhang; overhang
  polarity of EcoRI/BamHI/HindIII vs. PstI-type 3′ ends; origin-spanning
  EcoRI site). Compared with Biopython 1.88 `Bio.Restriction` (REBASE-
  derived data) on pUC19, NCBI L09137.2: enzyme definitions 22/22, digest
  cut positions and fragments 132/132, circular restriction map 132/132,
  end-compatibility matrix 484/484 (`docs/PHASE6_INITIAL_AUDIT.md` §3).
  Type IIS (2026-11-11, `tests/test_external_reference.py`): definitions
  4/4 and 32 digests of a 3,000 bp sequence with sites in both orientations
  and at both ends (linear and circular, four origin rotations), 140 cut
  positions, 0 differences. Not compared against NEBcutter.

## Primer design (`analyses/primer.py`)

Thermodynamic model (default):

- Nearest-neighbor ΔH/ΔS: SantaLucia (1998) unified parameters with
  terminal AT/GC initiation.
- Salt: ΔS correction 0.368·(N−1)·ln[Na⁺eq] (SantaLucia 1998).
- Divalent cations: [Na⁺eq] = [Na⁺] + 120·√([Mg²⁺] − [dNTP]) in mM
  (von Ahsen et al. 2001).
- Tm = ΔH / (ΔS + R·ln(C/4)) for non-self-complementary oligos
  (C/1 and −1.4 e.u. symmetry term for self-complementary ones).
- Conditions: 50 mM Na⁺, 1.5 mM Mg²⁺, 0.2 mM dNTP, 250 nM primer.

Alternatives: `nn_parameter_set` = `santalucia_2004` or `synapse_legacy`;
`salt_model` = `na_only` or `synapse_legacy_4x` (the previous rule
[Na⁺] + 4·([Mg²⁺] − [dNTP])).

Design defaults: primer 18–25 nt (optimum 20), GC 40–60%, Tm 57–63 °C
(optimum 60), ΔTm ≤ 2 °C, amplicon 80–1000 bp, homopolymer ≤ 4.

Coordinates: the menu takes the target region as 1-based inclusive
positions and prints primer and target positions the same way. The
design functions, `core.api` and exported results keep 0-based half-open
`start`/`end` (`sequence[start:end]` is the primer, reverse-complemented
for the reverse primer).

Validation:

- Known-answer test: the worked example in the Primer3 manual
  (CGTGACGTGACGGACT, 50 mM Na⁺, 50 nM) – ΔH and ΔS exact, Tm 50.550 °C vs.
  50.554 °C (the manual uses R = 1.987).
- Implementation cross-check against an independent implementation of the
  published equations (6 primers × 4 buffer conditions, agreement < 1e‑6 °C).
- External reference comparison (2026-11-02): primer3-py 2.4.0, which
  wraps the Primer3 C library and shares no code with Synapse, with
  `tm_method="santalucia"` and `salt_corrections_method="santalucia"`.
  81 primers (the manual example and 18-25-nt windows of a fixed test
  sequence) x 36 conditions (Na⁺ 50 and 100 mM; Mg²⁺ 0, 1.5 and 3 mM;
  dNTP 0, 0.2 and 0.8 mM; oligo 50 and 250 nM): 2,916 comparisons,
  max |ΔTm| 0.0032 °C, mean 0.0024 °C. The systematic offset comes from
  the gas constant (1.98720 vs 1.987). A reduced grid runs in the suite
  (`Primer3ReferenceComparison`, skipped without the `reference` extra).
- Hairpin and dimer screens are heuristics: the length of the longest
  perfectly complementary stem (≥ 4 / 7 / 10 bp: MODERATE / HIGH / VERY
  HIGH) and the 3′ complementarity, not free energies. Compared on
  2026-11-11 with primer3-py 2.4.0 (`calc_hairpin`, `calc_homodimer`,
  `calc_heterodimer`, same buffer) on 250 primers and 100 pairs
  (`tests/test_primer_structure_reference.py`): the stem length orders
  Primer3's ΔG for hairpins (Spearman ρ −0.78; median ΔG 0.00 / −1.77 /
  −4.82 kcal/mol for LOW / MODERATE / HIGH) and self-dimers (ρ −0.74), and
  less well for heterodimers (ρ −0.51). The classes overlap: a MODERATE
  self-dimer can reach −15.6 kcal/mol, and VERY HIGH (driven by 3′
  complementarity) is not more stable than HIGH in the median. A severity
  is therefore a screening flag; Primer3's ΔG is the better measure of
  structure stability.
- Predicted Tm is not an annealing temperature, and agreement with Primer3
  shows that the equations are implemented as published, not that the model
  predicts experiments.

## Mutation analysis (`analyses/mutation.py`)

- Global alignment (Needleman–Wunsch; match 2, mismatch −1, gap −2); pure
  length changes use the common prefix/suffix representation.
- Coding consequences assume that position 1 is the first base of codon 1.
  A CDS can be taken from the selected ORF; plausibility warnings are given
  when that assumption looks wrong. Start-loss is reported only for codon 1.
- Notation is HGVS-style (den Dunnen et al. 2016), not checked against an
  HGVS validator. `c.` positions count bases of the supplied sequence and
  equal HGVS c. numbering only if that sequence is the CDS. Protein changes
  are predicted from the DNA and therefore written in parentheses, with
  one-letter codes: `p.(R12C)` missense, `p.(R12=)` synonymous, `p.(W24*)`
  nonsense, `p.(M1?)` for a changed initiator codon (effect on initiation
  unknown), `p.(*110Qext*?)` for a lost stop whose new stop is not
  determined. The mutation report shows the notation next to the effect,
  and the variant bridge stores it as `protein_change` (2026-11-11).
- Comparisons above 2,000 nt are flagged (O(n·m) time and memory).

## Variant analysis (`analyses/variant_analysis.py`)

- VCF 4.x text or gzip/BGZF. Multi-allelic records are split; INFO/FORMAT
  values follow the header `Number=` (A, R, G; VCF 4.3 reserved keys as
  fallback). Number=G uses the VCF genotype order F(j,k) = k(k+1)/2 + j.
  Other ALT alleles in GT become 0, as in `bcftools norm -m-`.
- Normalization: trimming and left-shifting through repeats (Tan et al.
  2015) when the contig sequence is supplied; REF alleles are checked
  against it (`REF_MISMATCH`).
- Contig names: CHROM is matched to the reference names exactly by default
  (`contig_matching="exact"`). `contig_matching="chr_prefix"` also pairs
  `1` with `chr1`; use it only when both describe the same assembly.
  Mitochondrial names are never paired (hg19 chrM, 16,571 bp, is not the
  rCRS MT of GRCh37/GRCh38, 16,569 bp), case variants are not folded, and
  CHROM is never renamed. A variant whose contig has no reference is only
  trimmed and is flagged `REFERENCE_CONTIG_NOT_FOUND`; an aliased one is
  flagged `CONTIG_NAME_ALIASED`. Contig lengths are not compared, and
  RefSeq/GenBank accession names are not translated.
- FILTER `.` counts as "not filtered", not as failed.
- Validation: known-answer tests; the left-alignment algorithm was checked
  by hand against Tan et al. (2015) on constructed repeat cases
  (`docs/PHASE6_INITIAL_AUDIT.md` §3.7); not compared against
  `bcftools norm` (no native Windows build; recommended next validation
  step).

## CRISPR guide design (`analyses/crispr.py`)

- PAM search on both strands (SpCas9 NGG, SpCas9-NG, SaCas9 NNGRRT,
  Cas12a TTTV, …). Restricting to an ORF uses the ORF's genomic region.
- Ranking is a transparent weighted sequence heuristic, not a measured or
  predicted editing efficiency.
- Off-target search covers only sequences supplied to Synapse, not a genome.
- Recognition-site scanning is cached on (sequence, recognition sequence,
  topology), the three values it depends on. This removes the cloning
  planner's repeated rescanning (S18) and changes no result.
- Cut site (D8): SpCas9/SpCas9-NG blunt 3 nt 5' of the PAM (after guide
  position 17; Jinek et al. 2012), SaCas9 blunt 3 nt 5' of the PAM
  (after position 18; modeled as SpCas9, Ran et al. 2015), Cas12a
  staggered after position 18 (PAM strand) and 23 (target strand), 5-nt
  5' overhang (Zetsche et al. 2015). Reported as forward-strand
  positions: the break lies between p and p+1. Geometry only, not
  cleavage efficiency.
- The guide's own site in a supplied local reference is reported as
  "on-target" and not counted as an off-target (D9).
- Guide seed (S23): the seed is the PAM-proximal end of the guide in
  targeting orientation, i.e. the 3′ end for 3′-PAM systems (Cas9) and the
  5′ end for 5′-PAM systems (Cas12a). The PAM-proximal seed concept is
  established for Cas9 (Sternberg et al. 2014, PMID 24476820); the default
  length of 8 nt is the conservative end of the commonly used 8–12 nt range
  and is a Synapse convention, not a value prescribed by that work.
  `pam_proximal_seed()` is the single source of this rule for both the
  stored seed and the local specificity screen.
- Screening thresholds (S4): the GC range (30–70 %), the homopolymer limit
  (4) and the local-mismatch limit (3) are Synapse screening conventions
  consistent with common guide-design practice. They are not published
  cut-offs, they are configurable through `DesignConfig`, and they are
  quality filters rather than efficiency predictions.
- Poly-T warning (S5): a run of ≥ 4 T in a guide is reported as a *possible*
  RNA polymerase III termination signal, which is relevant only when the
  guide is expressed from a U6 or H1 promoter; synthetic or in-vitro
  transcribed guides are unaffected. The threshold is deliberately
  conservative. Verosloff et al. 2021 (PMID 33811918) found four uracils
  sufficient to terminate at Pol III type 2 promoters, whereas type 3
  promoters such as U6 required longer tracts, or shorter tracts combined
  with upstream RNA secondary structure. Synapse therefore screens for the
  motif and does not predict that termination occurs.
- Base-editing windows (S4): protospacer positions, 1-based from the
  PAM-distal end, inclusive. CBE (4, 8) follows the window reported for BE3
  (Komor et al. 2016, PMID 27096365); ABE (4, 8) is a conservative superset
  of the narrower window reported for ABE7.10 (Gaudelli et al. 2017,
  PMID 29160308). The `CBE+`/`ABE+` variants widen the window by one
  position and are Synapse screening conventions. A position inside the
  window is a candidate for editing, not an expected edit.
- Prime editing (S3): PBS and RTT are sequence-level constructions. The RTT
  is the **reference** template spanning the edit site; it does not encode
  the alternative allele and is therefore not a ready-to-use pegRNA RTT.
  The reported `edit_in_rtt` states that the RTT window covers the edit
  position — previously it only compared the RTT length with the requested
  length, which said nothing about the edit (prime editing after Anzalone
  et al. 2019, PMID 31634902). No thermodynamic, structural or efficacy
  model is applied.

## BLAST (`analyses/blast.py`)

- Remote NCBI BLAST (after an explicit privacy confirmation) or local BLAST+.
- Query coverage = union of HSP query intervals / query length (NCBI
  definition); the per-HSP value is reported separately.
- Evidence bands are assigned only to hits with E ≤ 1e‑3
  (`DEFAULT_EVIDENCE_EVALUE_THRESHOLD`, configurable); the threshold is a
  screening convention, not a universal significance level.
- Stored: parameters reported by BLAST, database title/date/version.
- A remote search can be collected later from its request ID (P1-7).
  Only the RID is sent; BLAST output does not contain the query
  sequence, so a resumed run keeps the reported query length and an
  empty sequence rather than an invented one, and is marked as resumed.
- Query type: from `core.alphabets` (declared, or inferred from the
  letters). An undetermined query is not guessed: the menu asks for DNA
  or protein, the programmatic entry points refuse it. RNA is searched as
  DNA with U read as T, recorded in the query's notes. The protein
  alphabet includes J, U and O, which BLAST+ accepts.
- Tested: parsing and scoring with synthetic XML; local blastp against
  small protein databases built with `makeblastdb` when BLAST+ is
  installed. Not tested: live remote searches and large local databases.

## Literature (`analyses/literature.py`)

- PubMed, Europe PMC, OpenAlex, Crossref; each source may fail independently.
- Requests are throttled per host (NCBI 3/s, 10/s with an API key) and
  retried on HTTP 429/5xx (max. 3 attempts, `Retry-After` honored).
- Relevance scores measure text matching only, not evidence quality.

## Reproducibility

- Every stored result is an `AnalysisResult` (parameters, input SHA-256,
  warnings, limitations, timestamp, Synapse and Python version).
- Menu `[R]` exports all stored results (JSON complete, Markdown summary).
- `core/api.py` offers the same analyses without interactive input.

## References

- Cornish-Bowden A (1985) Nomenclature for incompletely specified bases in
  nucleic acid sequences: recommendations 1984. Nucleic Acids Res
  13:3021–3030. PMID 2582368.
- IUPAC-IUB Joint Commission on Biochemical Nomenclature (1984)
  Nomenclature and symbolism for amino acids and peptides. Eur J Biochem
  138:9–37. PMID 6692818.
- Jinek M et al. (2012) A programmable dual-RNA-guided DNA endonuclease
  in adaptive bacterial immunity. Science 337:816–821. PMID 22745249.
- Ran FA et al. (2015) In vivo genome editing using Staphylococcus aureus
  Cas9. Nature 520:186–191. PMID 25830891.
- Zetsche B et al. (2015) Cpf1 is a single RNA-guided endonuclease of a
  class 2 CRISPR-Cas system. Cell 163:759–771. PMID 26422227.
- Nishimasu H et al. (2018) Engineered CRISPR-Cas9 nuclease with expanded
  targeting space. Science 361:1259–1262. PMID 30166441.
- Sternberg SH et al. (2014) DNA interrogation by the CRISPR RNA-guided
  endonuclease Cas9. Nature 507:62–67. PMID 24476820.
- Komor AC et al. (2016) Programmable editing of a target base in genomic
  DNA without double-stranded DNA cleavage. Nature 533:420–424.
  PMID 27096365.
- Gaudelli NM et al. (2017) Programmable base editing of A•T to G•C in
  genomic DNA without DNA cleavage. Nature 551:464–471. PMID 29160308.
- Anzalone AV et al. (2019) Search-and-replace genome editing without
  double-strand breaks or donor DNA. Nature 576:149–157. PMID 31634902.
- Verosloff MS et al. (2021) RNA sequence and structure determinants of
  Pol III transcriptional termination in human cells. J Mol Biol
  433:166978. PMID 33811918.
- SantaLucia J Jr (1998) A unified view of polymer, dumbbell, and
  oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95:1460–1465.
- SantaLucia J Jr, Hicks D (2004) The thermodynamics of DNA structural
  motifs. Annu Rev Biophys Biomol Struct 33:415–440.
- von Ahsen N, Wittwer CT, Schütz E (2001) Oligonucleotide melting
  temperatures under PCR conditions. Clin Chem 47:1956–1961.
- Tan A, Abecasis GR, Kang HM (2015) Unified representation of genetic
  variants. Bioinformatics 31:2202–2204.
- The Variant Call Format Specification, VCFv4.3 (samtools/hts-specs).
- den Dunnen JT et al. (2016) HGVS recommendations for the description of
  sequence variants: 2016 update. Hum Mutat 37:564–569.
- Needleman SB, Wunsch CD (1970) A general method applicable to the search
  for similarities in the amino acid sequence of two proteins. J Mol Biol
  48:443–453.
- Camacho C et al. (2009) BLAST+: architecture and applications. BMC
  Bioinformatics 10:421.
- NCBI: The Genetic Codes; ORFfinder documentation; Primer3 manual
  (`PRIMER_TM_FORMULA`, `PRIMER_SALT_CORRECTIONS`).
- Gardiner-Garden M, Frommer M (1987) CpG islands in vertebrate genomes.
  J Mol Biol 196:261–282. (Bibliographic details not verified at the
  source in Phase 6.)
- Guruprasad K, Reddy BVB, Pandit MW (1990) Correlation between stability
  of a protein and its dipeptide composition: a novel approach for
  predicting in vivo stability of a protein from its primary sequence.
  Protein Eng 4:155–161. (Journal, volume and pages as cited in
  Biopython's `ProtParamData`; not verified at the source in Phase 6.)
- Cock PJA et al. (2009) Biopython: freely available Python tools for
  computational molecular biology and bioinformatics. Bioinformatics
  25:1422–1423. (Reference implementation in the tests, version 1.88.
  Bibliographic details not verified at the source in Phase 6.)
- Untergasser A et al. (2012) Primer3 - new capabilities and interfaces.
  Nucleic Acids Res 40:e115. (primer3-py 2.4.0 wraps Primer3 in the tests.
  Bibliographic details not verified at the source in Phase 6.)
