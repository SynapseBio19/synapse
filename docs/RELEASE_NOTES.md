# Synapse release notes

## 0.3.0 (prepared, not yet released)

The first version meant for users outside the project. The code in the
repository is complete for this release; the release itself waits for the
owner's decisions listed under "Before the first release". The detailed,
dated record of every change is `docs/CHANGES.md`.

### What Synapse offers

- **Twelve analyses in one session**: sequence overview, ORF analysis,
  RNA analysis, protein analysis, restriction and cloning, primer design,
  mutation and variant analysis, CRISPR guide design (partial), BLAST
  (NCBI or local BLAST+), literature search (PubMed, Europe PMC,
  OpenAlex, Crossref), session report. The AI Experiment Planner is shown
  as Coming Soon.
- **All 27 NCBI genetic codes**, carried consistently from the ORF to the
  RNA codon view, the protein and the report.
- **HGVS notation** for DNA and predicted protein changes (`c.215C>G`,
  `p.(P72R)`, 3' rule, dup, delins).
- **Complete bacterial genomes**: every analysis runs on a 4.6 Mb genome
  (docs/BENCHMARKS.md).
- **Traceable results**: every result records the input checksum, all
  parameters, the genetic code, database and enzyme-table versions,
  warnings and limitations; JSON and Markdown reports.
- **Python interface** (`core.api`) for every analysis, without the menus.
- **`synapse` command** after `pip install .`; the core needs only the
  Python standard library (3.10 or newer).
- **User guide** for students and laboratory users (docs/USER_GUIDE.md).

### Validation

Compared with: NCBI ORFfinder (69/69 ORFs), the RefSeq annotation of
*E. coli* K-12 (every annotated gene of at least 90 nt found at its stop
codon), NCBI's genetic-code tables (27/27), Primer3 (Tm within
0.003 °C), REBASE and Biopython (restriction and Type IIS cuts),
bcftools norm (29/29), Mutalyzer (11/11) and 1,500 ClinVar variants
(c.: 1,489 identical and 11 equivalent; p.: 1,186/1,186), ExPASy ProtParam
(399/400 instability weights), ViennaRNA, and NCBI BLAST with the
refseq_rna database. Synapse is not experimentally or clinically
validated (docs/METHODS.md, docs/BENCHMARKS.md).

### Defects found and fixed during the preparation

Selected from docs/CHANGES.md; each fix has a test that failed before it.

- BLAST: an E-value of 0 (a perfect hit) was read as infinity, so the
  best hits counted as not significant and were sorted last (local and
  NCBI searches).
- BLAST: temporary database views were never removed, and "Remove
  Database" left one file behind.
- HGVS: DNA changes used a non-standard notation; in an allele with
  several changes, the 3' rule could move one change onto another
  (found with the ClinVar benchmark).
- Variants: VCF records without an alternative allele were dropped
  silently; they are now reported.
- Sequence overview: one statistic depended on Python's hash seed, so the
  same input could give slightly different values in two runs.
- Reports: RNA and protein results did not state the genetic code.
- Literature: BibTeX exports broke LaTeX on titles with `&` or `%`.
- Installation: an installed Synapse would have written logs, exports and
  databases into the Python installation.
- Genome-length input: the ORF search took 8 minutes and a full-panel
  restriction digest did not finish in 15 minutes on *E. coli*; the
  overview needed 3.5 GB of memory. Now about 10 s, 25 s and well under
  1 GB, with identical results.
- Terminal: three protein errors lacked the ERROR prefix; Back inside a
  cloning analysis left the cloning menu.

### Known limitations

- CRISPR: no CFD off-target score and no genome-wide off-target search;
  CRISPR development is paused by owner decision.
- The global alignment of the mutation module uses linear gap costs, so a
  long deletion-insertion can be written as several changes (an
  equivalent but not the preferred HGVS description).
- The hairpin screen has no energy model; the optional ViennaRNA view does.
- Hairpin and primer-dimer screens are heuristics, not free energies.
- Tested on Windows with Python 3.10 and 3.14. The CI configuration adds
  Linux with Python 3.10-3.14 but has not run on GitHub yet; macOS is not
  tested.

### Before the first release (owner)

1. ~~Choose the licence.~~ Done 2026-11-20: AGPL-3.0-or-later, copyright
   Ali Onur Bayraktar (`LICENSE`, `pyproject.toml`, `CITATION.cff`,
   README).
2. Create the GitHub repository and push `main`. The workflow
   `.github/workflows/ci.yml` then runs the tests on Linux (Python
   3.10-3.14) and Windows. Fix what the first run finds before tagging.
3. Set the version to `0.3.0` in `pyproject.toml`, `core/models.py`
   (`SYNAPSE_VERSION`) and `CITATION.cff`; add `date-released` and
   `repository-code` to `CITATION.cff`.
4. Tag the commit `v0.3.0` and create a GitHub release with these notes.
5. Optional: archive the release on Zenodo for a DOI; publish on PyPI as
   `synapse-gen`.
