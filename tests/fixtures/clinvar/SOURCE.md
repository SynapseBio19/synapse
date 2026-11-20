# ClinVar HGVS benchmark set

`hgvs_benchmark.json` holds 1,500 ClinVar variant names with their RefSeq
coding sequences. It is the reference data of
`tests/test_hgvs_clinvar_benchmark.py`, which checks Synapse's HGVS c. and
p. descriptions on curated variants.

- Built by `scripts/benchmarks/clinvar_hgvs.py` on 2026-11-20 from NCBI
  E-utilities: ClinVar `esearch` (`<gene>[gene]`) and `esummary` (the
  variant title, e.g. `NM_000546.6(TP53):c.215C>G (p.Pro72Arg)`), and
  nuccore `efetch` (`rettype=fasta_cds_na`) for the coding sequences.
- Selection: for each gene the transcript named most often in ClinVar;
  exonic coding changes only (no intronic offsets, no UTR positions), the
  250 with the lowest ClinVar UIDs per gene. No clinical significance or
  submitter data is stored.

| Gene | Transcript | CDS (nt) | Variants |
|---|---|---:|---:|
| HBB | NM_000518.5 | 444 | 250 |
| KRAS | NM_004985.5 | 567 | 250 |
| VHL | NM_000551.4 | 642 | 250 |
| TP53 | NM_000546.6 | 1,182 | 250 |
| PTEN | NM_000314.8 | 1,212 | 250 |
| MECP2 | NM_001110792.2 | 1,497 | 250 |

- Licence: ClinVar and RefSeq records are produced by NCBI, an agency of
  the U.S. government, and are in the public domain
  (<https://www.ncbi.nlm.nih.gov/home/about/policies/>).
- The set tests notation only. Synapse does not interpret pathogenicity.

SHA-256 (checked by the test):

```
0348314292321b9646ff921505625a92b461e724718f53b18578d74e8df4a0c4  hgvs_benchmark.json
```
