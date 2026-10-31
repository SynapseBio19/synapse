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
  type and confidence into provenance and the report. Without a declared
  type the letters decide (see "Sequence types").
- Main menu: a record whose type cannot be determined is not guessed; the
  user declares DNA, RNA or protein. RNA records are analyzed in the DNA
  alphabet (U read as T), which the report records. Protein records open
  the protein menu (protein analysis, BLAST, report).
- `core.api`: nucleotide analyses accept DNA and RNA (RNA converted, with a
  warning); restriction digests and CRISPR guide design accept DNA only,
  because both act on DNA; `analyze_protein` takes proteins. Every entry
  point accepts `sequence_type=` for an undetermined sequence.
- Provenance: file SHA-256, sequence SHA-256, path, format, record index.
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
- Classification without a declared type (`classify_sequence`):
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
- The "heuristic score" is a ranking aid, not a probability of coding.
- Not implemented: 5′-partial ORFs without a start codon, ORFs across the
  origin of circular sequences.
- Validation: known-answer tests for frames, strands, coordinates and the
  table-2 read-through of TGA; not compared against ORFfinder output.

## RNA and protein (`analyses/rna.py`, `analyses/protein.py`)

- Translation in frame +1; ambiguous codons give `?` unless all expansions
  encode one amino acid. With `initiator_start=True` (protein analysis of
  a selected ORF, `core.api.analyze_protein_from_cds`) the first codon is
  decoded as Met only if it is an initiation codon of the genetic code —
  the same rule as the ORF translator.
- Protein analysis uses the genetic code of the selected ORF. The RNA
  module's three-frame codon view always uses table 1.
- Physicochemical values (mass, pI, GRAVY and similar) are sequence-based
  estimates. Checked against ExPASy ProtParam for ubiquitin (MW 8564.84 Da,
  pI 6.56, GRAVY −0.489, ε 1490 M⁻¹cm⁻¹); otherwise not independently
  validated.
- Unresolved residues ("?", from ambiguous codons) keep their position.
  Sliding windows (hydropathy, flexibility, transmembrane screen, signal
  peptide) and the instability index skip the windows or dipeptides that
  contain one instead of deleting the residue, so reported positions always
  refer to the real protein. Composition-based values (pI, GRAVY, extinction
  coefficient and similar) are computed over the resolved residues only.
- The transmembrane and signal-peptide outputs are screens based on the
  Kyte–Doolittle scale, not predictions from TMHMM/Phobius or SignalP.
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
- End compatibility requires the same overhang sequence *and* polarity
  (5′ vs 3′); blunt ends are compatible with blunt ends.
- Methylation sensitivity is shown as a note (ClaI, XbaI, DpnI, HpaII);
  blocking by methylation and star activity are not modelled.
- Validation: known-answer tests (ClaI cut positions and overhang; overhang
  polarity of EcoRI/BamHI/HindIII vs. PstI-type 3′ ends; origin-spanning
  EcoRI site). Compared with Biopython 1.88 `Bio.Restriction` (REBASE-
  derived data) on pUC19, NCBI L09137.2: enzyme definitions 22/22, digest
  cut positions and fragments 132/132, circular restriction map 132/132,
  end-compatibility matrix 484/484 (`docs/PHASE6_INITIAL_AUDIT.md` §3).
  Not compared against NEBcutter.

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

Design defaults: primer 18–25 nt (optimum 20), GC 40–60 %, Tm 57–63 °C
(optimum 60), ΔTm ≤ 2 °C, amplicon 80–1000 bp, homopolymer ≤ 4.

Validation:

- Known-answer test: the worked example in the Primer3 manual
  (CGTGACGTGACGGACT, 50 mM Na⁺, 50 nM) – ΔH and ΔS exact, Tm 50.550 °C vs.
  50.554 °C (the manual uses R = 1.987).
- Implementation cross-check against an independent implementation of the
  published equations (6 primers × 4 buffer conditions, agreement < 1e‑6 °C).
- Not run: Primer3 or Biopython executables (not installable in the
  development environment). Hairpin/dimer screens are heuristic and not
  independently validated. Predicted Tm is not an annealing temperature.

## Mutation analysis (`analyses/mutation.py`)

- Global alignment (Needleman–Wunsch; match 2, mismatch −1, gap −2); pure
  length changes use the common prefix/suffix representation.
- Coding consequences assume that position 1 is the first base of codon 1.
  A CDS can be taken from the selected ORF; plausibility warnings are given
  when that assumption looks wrong. Start-loss is reported only for codon 1.
- Notation is HGVS-like (`c.`/`p.` prefixes on Synapse notation), not
  validated against the HGVS recommendations (den Dunnen et al. 2016).
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

## BLAST (`analyses/blast.py`)

- Remote NCBI BLAST (after an explicit privacy confirmation) or local BLAST+.
- Query coverage = union of HSP query intervals / query length (NCBI
  definition); the per-HSP value is reported separately.
- Evidence bands are assigned only to hits with E ≤ 1e‑3
  (`DEFAULT_EVIDENCE_EVALUE_THRESHOLD`, configurable); the threshold is a
  screening convention, not a universal significance level.
- Stored: parameters reported by BLAST, database title/date/version.
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
