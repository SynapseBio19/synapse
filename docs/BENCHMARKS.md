# Synapse benchmarks

Two kinds of measurement, both repeatable with the scripts named below:

1. **Genome scale**: does every analysis run on a complete bacterial
   genome, and how long does it take and how much memory does it need?
2. **Curated data**: do Synapse's results agree with curated reference
   annotations, not only with other programs?

All values are from 2026-11-20 unless stated otherwise.

## Machine and data

| Item | Value |
|---|---|
| CPU | Intel Core Ultra 7 255HX (20 logical processors) |
| Memory | 15.4 GB |
| System | Windows 11 (10.0.26200), CPython 3.14.6 |
| Genome | *Escherichia coli* K-12 substr. MG1655, RefSeq NC_000913.3, 4,641,652 bp, GC 50.79 % |
| Annotation | RefSeq GFF3 of NC_000913.3, 4,340 CDS |
| Variants | 1,500 ClinVar names on six RefSeq transcripts (`tests/fixtures/clinvar/`) |

Sources, licences and checksums: `docs/TOOL_AND_DATA_LEDGER.md`, section
"Phase 10". The genome and the annotation are downloaded once into
`.cache/benchmarks/` (not in Git) and checked by SHA-256.

## 1. Genome scale

`python scripts/benchmarks/genome_scale.py` runs each analysis on the
whole genome in its own process with a time limit, and reports wall time,
CPU time and the peak working set of that process. Wall times vary by a
few percent between runs.

| Analysis | Setting | Before (code of `a2898aa`) | After (Phase 10) | Peak memory after |
|---|---|---:|---:|---:|
| Load FASTA | 4.6 Mb record | 0.3 s | 0.3 s | 42 MB |
| Sequence overview | defaults (100-nt profiles, step 1) | 154 s | 41 s | 0.64 GB (3.5 GB before) |
| ORF search | defaults (ATG, ≥ 75 nt, longest per stop) | 499 s | 11.7 s | 0.50 GB |
| ORF search | code 11, ATG/GTG/TTG, ≥ 90 nt, 66,712 ORFs | not measured | 87 s | – |
| RNA analysis | whole genome as one transcript | > 600 s (time limit) | 137 s | 2.1 GB |
| Restriction digest | all 22 enzymes, 119,067 cuts, 93,240 fragments | > 900 s (time limit) | 5.9 s | 0.48 GB |
| Primer design | target 2,000,001-2,000,600 | 1.6 s | 1.7 s | 65 MB |
| Protein analysis | thrA (2,463 nt) | 0.35 s | 0.27 s | 42 MB |
| Sequence comparison | 2 Mb against 2 Mb with one SNP | 84 s | 73 s (unchanged code; run-to-run variation) | 0.45 GB |

**Identical results.** The faster code returns exactly what the former
code returned. Evidence:

- Differential tests compare every changed function with its former
  implementation on random input with IUPAC codes, lower case and edge
  sizes, float for float and in order: `tests/test_window_scans.py`,
  `test_rna_window_scans.py`, `test_orf_relationship_sweep.py`,
  `test_orf_feature_speed.py`, `test_restriction_fragment_index.py`.
- End to end, the code of `a2898aa` and the current code give
  bit-identical results (SHA-256 of the complete output) on *E. coli*
  segments:

| Case | Before | After | Result |
|---|---:|---:|---|
| Sequence overview, 150 kb | 3.1 s | 0.9 s | identical |
| ORF search, defaults, 400 kb | 7.7 s | 0.6 s | identical |
| ORF search, ATG/GTG/TTG, nested, ≥ 60 nt, 150 kb | 155.6 s | 20.5 s | identical |
| ORF search, circular, 60 kb | 3.0 s | 1.0 s | identical |
| RNA analysis, 60 kb | 11.3 s | 2.2 s | identical |
| Digest, 22 enzymes, linear, 400 kb | 18.2 s | 0.8 s | identical |
| Digest, 22 + 4 Type IIS enzymes, circular, 150 kb | 7.0 s | 0.4 s | identical |

**One documented behavior change, for long input only.** Above 200,000
profile windows (about 200 kb) the sequence overview and the RNA analysis
store every k-th window of their GC and complexity profiles and record
this in `profile_resolution`. The DNA profile statistics (minimum,
maximum, average and their positions) are still computed over every
window. Without the thinning, the two profiles of the genome took 3.5 GB.

**What still costs.**
- The RNA analysis of a whole genome keeps every hairpin-screen
  candidate, which accounts for most of its 2.1 GB. The RNA module is
  meant for transcripts; for a genome, analyze the ORFs.
- The global alignment of the mutation module is exact. Similar sequences
  run in a band around the diagonal, so 2 Mb with one difference takes
  about 1.5 minutes; unrelated long sequences take time proportional to
  the product of their lengths. The module recommends at most 2,000 nt.
- An ORF search with alternative start codons evaluates every in-frame
  start (about 110,000 candidates per Mb); 87 s for the genome.

## 2. Curated data

### ORFs against the RefSeq gene annotation

`python scripts/benchmarks/orf_refseq.py`

- **Truth set.** The 4,297 single-segment CDS of NC_000913.3 that are not
  pseudogenes and have no ribosomal-slippage exception (excluded: 18
  pseudogene and 3 slippage CDS); 3 of them are selenoprotein genes.
- **Synapse settings.** Code 11, start codons ATG, GTG and TTG, one ORF
  per stop codon (the longest), complete ORFs only, minimum 90 nt.
- **Criterion.** A CDS is found when an ORF on the same strand ends at the
  same base. The stop is unambiguous; the start is a choice among in-frame
  start codons, so it is evaluated separately.

| Result | Count |
|---|---:|
| Regular CDS found at their stop codon | 4,250 of 4,294 (98.98 %) |
| Missed, all shorter than the 90-nt minimum | 44 |
| CDS of at least 90 nt found | 4,250 of 4,250 (100 %) |
| ORF start identical to the annotated start | 1,470 (34.6 % of found) |
| ORF start upstream (a longer ORF) | 2,780 (65.4 %) |
| ORF start downstream of the annotated start | 0 |
| Selenoprotein genes | 3: the ORF ends at the internal UGA, as a plain ORF search must |

The longest-ORF rule picks the first in-frame start codon, so two thirds
of the ORFs begin upstream of the annotated start; no ORF begins inside
an annotated gene. Precision is not reported: an ORF search is not a gene
predictor, and long ORFs on the antisense strand of real genes are not
genes. Of the 9,135 ORFs of at least 300 nt, 3,846 end at an annotated
stop codon.

### HGVS notation against ClinVar

`python -m unittest tests.test_hgvs_clinvar_benchmark` (all 1,500 with
`SYNAPSE_CLINVAR_FULL=1`; the default run uses every 10th variant and the
documented cases).

| Gene | Transcript | CDS | Variants |
|---|---|---:|---:|
| HBB | NM_000518.5 | 444 | 250 |
| KRAS | NM_004985.5 | 567 | 250 |
| VHL | NM_000551.4 | 642 | 250 |
| TP53 | NM_000546.6 | 1,182 | 250 |
| PTEN | NM_000314.8 | 1,212 | 250 |
| MECP2 | NM_001110792.2 | 1,497 | 250 |

Change types: 1,173 substitutions, 223 deletions, 63 duplications, 12
insertions, 29 deletion-insertions. For each variant the observed coding
sequence is rebuilt from ClinVar's c. name, and Synapse compares it with
the reference (alignment, mutation calling, HGVS formatting).

| Result | Count |
|---|---:|
| c. description identical to ClinVar | 1,489 of 1,500 |
| c. description equivalent, written as several changes | 11 (deletion-insertions split by the alignment) |
| Descriptions that rebuild the observed sequence | 1,500 of 1,500 |
| Protein change of codon substitutions identical | 1,186 of 1,186 |

- **Defect found.** Before the fix, PTEN c.528_534delinsA was written as
  `c.[532_534del;531_532del;534del]`, three overlapping deletions that
  describe no sequence (docs/CHANGES.md, Phase 10 B).
- **Notation policy.** ClinVar writes `p.Met1Arg` and `p.Ter214Leu`;
  Synapse writes the HGVS recommendations `p.(M1?)` (effect on initiation
  unknown) and `p.(*214Lext*?)` (extension, new stop undetermined). The
  benchmark maps one to the other.
- **The 11 split deletion-insertions.** The alignment uses linear gap
  costs, so a long deletion-insertion can come out as several small
  changes (for example `c.[217_218del;220_221del]` for ClinVar's
  `c.217_221delinsT`). Both describe the same sequence; ClinVar's is the
  preferred HGVS form. Affine gap costs would avoid this but change the
  alignment model, so they are not introduced without a separate decision.

### Not done

- **Primer melting temperatures against measured values.** No published
  dataset of measured Tm values with a clear licence was obtained, and
  values were not reconstructed from memory. Synapse's Tm agrees with
  Primer3 within 0.003 °C (`docs/METHODS.md`), which shows that the
  published model is implemented, not that the model predicts a
  particular laboratory's melting curve.

## 3. Live NCBI BLAST

One remote search with a public sequence (human GAPDH, NM_002046.7
nt 301-700) against `refseq_rna` on 2026-11-20.

- Submitted with `core.api.blast_search(..., engine="remote",
  send_to_ncbi=True)`: NCBI accepted it as RID FZ5HPRV0014.
- After the default 300 s wait Synapse stopped and reported the RID and
  how to collect it, as designed.
- `core.api.blast_resume` then waited another 2,400 s; NCBI still had
  not finished, and the resume ended with the same clear message.
- **No result was obtained**, as on 2026-11-19. The live confirmation of
  the remote path against `refseq_rna` is therefore still open. The remote
  path itself was confirmed live on 2026-11-13 (`nt`) and 2026-11-17
  (`core_nt`), and its result parser is the one checked on the local
  `refseq_rna` searches (2026-11-18). Retry when NCBI is less busy, for
  example on a weekend morning (BLAST menu, entry 2, with the RID while
  NCBI keeps it).

## Repeating the benchmarks

```
python scripts/benchmarks/genome_scale.py --json genome.json
python scripts/benchmarks/orf_refseq.py --json orf_refseq.json
SYNAPSE_CLINVAR_FULL=1 python -m unittest tests.test_hgvs_clinvar_benchmark
python scripts/benchmarks/clinvar_hgvs.py      # rebuilds the ClinVar set (network)
```

The first two download the genome and annotation once (about 7 MB).
`clinvar_hgvs.py` queries NCBI and writes a new fixture; the test checks
the fixture's checksum, so a rebuilt set is a deliberate change.
