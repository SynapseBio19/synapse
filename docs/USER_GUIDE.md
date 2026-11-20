# Synapse user guide

This guide is for students, teachers and laboratory scientists who want
to analyze DNA, RNA and protein sequences with Synapse. It explains the
concepts behind each analysis, walks through every menu with a real
example, and says how far each result can be trusted. No programming
knowledge is needed for the menus. The Python interface is described at
the end.

The technical reference for every method, parameter and validation is
[docs/METHODS.md](METHODS.md). This guide does not repeat it; it links to
it where a detail matters.

## Contents

1. [What Synapse is, and what it is not](#1-what-synapse-is-and-what-it-is-not)
2. [Installation](#2-installation)
3. [Starting a session and loading a sequence](#3-starting-a-session-and-loading-a-sequence)
4. [How to read the screens](#4-how-to-read-the-screens)
5. [Tutorial: the plasmid pUC19 from start to report](#5-tutorial-the-plasmid-puc19-from-start-to-report)
6. [Module reference](#6-module-reference)
7. [Reports, exports and reproducibility](#7-reports-exports-and-reproducibility)
8. [Using Synapse from Python](#8-using-synapse-from-python)
9. [How far the results can be trusted](#9-how-far-the-results-can-be-trusted)
10. [Troubleshooting](#10-troubleshooting)
11. [Glossary](#11-glossary)

---

## 1. What Synapse is, and what it is not

Synapse is a sequence-analysis workbench for the terminal. It brings
twelve analyses together in one session: sequence statistics, open
reading frames, RNA and protein properties, restriction digests and
cloning, primer design, mutation and variant analysis, CRISPR guide
design, BLAST and a literature search. A sequence loaded once flows
through all of them; the ORF you select in the ORF module is the one the
RNA and protein modules analyze.

Synapse is built for **research support and teaching**:

- Every result states the method, the parameters and its limitations.
- Heuristic scores are called heuristic; no score is presented as a
  probability.
- The calculations were compared with established tools and curated data
  (Primer3, NCBI ORFfinder, bcftools, Mutalyzer, ExPASy ProtParam,
  ViennaRNA, Biopython, NCBI BLAST, the RefSeq annotation of *E. coli*
  and 1,500 ClinVar variants; see section 9).

Synapse is **not**:

- a replacement for Benchling, Geneious, SnapGene or AlphaFold;
- a clinical tool: it does not assess the pathogenicity of variants;
- experimentally validated: whether a primer works in a PCR, or a guide
  RNA cuts, can only be shown in the laboratory.

Two parts are deliberately incomplete:

- **CRISPR Guide Design** finds guides and screens off-targets in the
  sequences you supply. It does not compute CFD scores and does not search
  a whole genome. Treat its output as a candidate list.
- **AI Experiment Planner** is shown as *Coming Soon*. It is planned but
  not available; selecting it runs nothing.

## 2. Installation

Synapse needs **Python 3.10 or newer** (tested on 3.10 and 3.14 on
Windows). The analyses themselves use only the Python standard library.

### From the project folder

Open a terminal (on Windows: PowerShell) in the Synapse folder and run:

```
python -m venv .venv
.venv\Scripts\activate          (Windows)
source .venv/bin/activate       (Linux, macOS)
pip install .
```

Then start Synapse with the command `synapse`, or, without installing,
with `python main.py` from the project folder.

### Optional components

| Install | Adds |
|---|---|
| `pip install ".[structure]"` | ViennaRNA: the energy-based RNA secondary-structure view |
| `pip install ".[reference]"` | Biopython and primer3-py: only for the comparison tests |
| BLAST+ from NCBI | Local BLAST searches (Synapse finds it on the PATH or via `SYNAPSE_BLAST_BIN`) |

On Windows 11, *Smart App Control* can block the compiled parts of
ViennaRNA and primer3-py. Synapse then says that the view is unavailable
and every other analysis keeps working.

### Where Synapse stores files

From the project folder, logs go to `logs/`, exports to `exports/` and
local BLAST databases to `data/blast_db/`. An installed Synapse uses a
folder in your user profile. `SYNAPSE_HOME` moves all of them, and
`SYNAPSE_BLAST_DB_ROOT` moves the BLAST databases alone. BLAST databases
are large (refseq_rna is about 81 GB): keep them outside folders that
OneDrive or another cloud service synchronizes.

## 3. Starting a session and loading a sequence

After the start Synapse asks for a sequence file:

```
======================================================================
LOAD SEQUENCE FILE
======================================================================
Formats:               FASTA (.fasta/.fa) or GenBank (.gb/.gbk)
Topology:              FASTA headers may contain topology=circular
Example:               ...\data\examples\dna.fasta
NOTE: Press Enter to load the example; type Q to quit.
```

- **Enter** loads the bundled example `data/examples/dna.fasta`: a 99-bp
  test gene and thirteen variants of it (substitutions, insertions,
  deletions, frameshifts), made for the mutation module.
- Type a path to load your own file. `data/examples/pUC19.fasta` holds
  the cloning vector pUC19 (GenBank L09137.2, 2,686 bp, circular), used in
  the tutorial below.
- **Q** quits.

**File formats.** FASTA files may hold several records. GenBank files
also carry the topology (linear or circular) and the molecule type. In a
FASTA header, `topology=circular` marks a circular molecule such as a
plasmid. Circularity matters: on a circle, restriction fragments and ORFs
can run across the origin.

**Sequence type.** Synapse decides whether a record is DNA, RNA or
protein from the file when it says so, and otherwise from the letters.
When the letters cannot decide (for example `ACGT...` could also be a
peptide), Synapse asks instead of guessing.

**Choosing a record.** The list shows each record with its length and
topology. Choose one by its number; **F** loads another file. The
sequence overview opens automatically, then the main menu:

```
SEQUENCE ANALYSIS
[1] Sequence Overview
[2] ORF Analysis
[3] RNA Analysis
[4] Protein Analysis

DESIGN AND COMPARISON
[5] Restriction & Cloning Analysis
[6] Primer Design
[7] Mutation & Variant Analysis
[8] CRISPR Guide Design

SEARCH AND PLANNING
[9] BLAST Search
[10] Literature Search
[11] AI Experiment Planner (Coming Soon)

SESSION
[P] Pipeline Status & Stored Results
[R] Export Analysis Report (JSON / Markdown)
[S] Select Another Sequence

[0] Exit
```

## 4. How to read the screens

- **Prefixes.** `ERROR:` means the step did not run (for example an
  unknown enzyme name); `WARNING:` means it ran but something deserves
  attention; `NOTE:` gives context. Synapse never fails silently.
- **Positions are 1-based and inclusive** on screen: position 1 is the
  first base, and 10-20 includes both ends. This is the convention of
  GenBank, NCBI and most biologists. Exported data and the Python
  interface use 0-based positions, as documented in METHODS.md.
- **Strands.** `+` is the sequence as loaded (5' to 3'); `-` is its
  reverse complement. Positions of `-` features are still given on the
  `+` strand coordinates.
- **0 goes back** one level; Enter at a pause continues.
- **Long output** is shortened with `...` in previews; the reports and
  exports hold the complete data.

## 5. Tutorial: the plasmid pUC19 from start to report

pUC19 is one of the most used cloning vectors. Its sequence carries an
ampicillin-resistance gene (*bla*), a fragment of *lacZ* with a multiple
cloning site, and an origin of replication. This tutorial finds them
with Synapse. Start Synapse, type `data/examples/pUC19.fasta` at the file
prompt and choose record 1.

### Step 1: Sequence overview

The overview opens by itself. For pUC19 it reports 2,686 bp, a GC content
of 50.63 %, the base composition, CpG statistics and sliding-window
profiles.

> **What is GC content?** The share of G and C among the bases. G-C
> pairs are held by three hydrogen bonds, A-T pairs by two, so GC-rich DNA
> melts at higher temperatures. Genomes differ widely (*E. coli* about
> 51 %, *Plasmodium falciparum* about 19 %), and GC-rich regions are
> harder to amplify by PCR.

> **What is a CpG island?** A stretch where the dinucleotide CG occurs
> about as often as its base composition predicts. In vertebrate genomes
> CG is depleted except in such islands, which often mark promoters.
> Synapse uses the classical screening criteria of Gardiner-Garden and
> Frommer (1987): at least 200 bp, GC at least 50 %, observed/expected CpG
> at least 0.6. They are a screen, not a proof of a regulatory element.

The **GC and complexity profiles** compute the value in a 100-bp window
moved base by base; the minimum and maximum tell you where unusual
regions are. For sequences with more than 200,000 windows (above about
200 kb), the stored profile keeps every k-th window to save memory, and
the overview says so; the statistics still use every window.

### Step 2: ORF analysis

Choose **2** in the main menu.

> **What is an ORF?** An open reading frame is a stretch of codons that
> starts with a start codon (usually ATG) and runs to the next in-frame
> stop codon (TAA, TAG or TGA in the standard code). DNA can be read in
> three frames on each strand, so six frames are searched. An ORF is a
> sequence feature: long ORFs are likely genes in bacteria, but an ORF
> alone does not prove that a protein is made.

The ORF menu:

```
[1] ORF Discovery      [2] ORF Candidates     [3] Compare ORFs
[4] ORF Relationships  [5] Select ORF         [6] Search Settings
```

Open **6 Search Settings** first. For a bacterial plasmid choose genetic
code **11** (bacterial, archaeal and plant plastid). Synapse knows all 27
NCBI genetic codes; the code changes which codons start and stop an ORF
and how codons translate. The defaults are a minimum length of 75 nt
(stop codon included), ATG as the only start codon, and one ORF per stop
codon (the longest). These defaults reproduce NCBI ORFfinder exactly,
except that ORFfinder's minimum length does not count the stop codon
(see METHODS.md).

**1 ORF Discovery** lists the ORFs. For pUC19 the longest one lies on the
minus strand from 1,626 to 2,486 (861 nt, 286 amino acids). That is
exactly the annotated *bla* gene of pUC19: Synapse found the
ampicillin-resistance gene without being told where it is.

**5 Select ORF** hands one ORF to the RNA and protein modules. Select the
861-nt ORF now.

### Step 3: RNA analysis

Choose **3**. Because an ORF is selected, the RNA module analyses its
transcript (without a selected ORF it would use the whole sequence; it
never switches silently). The menu offers the RNA overview, translation
and codon usage in the ORF's genetic code, RNA features (hairpin screen,
homopolymers, start/stop landscape), the RNA sequence and, with
ViennaRNA installed, **6 Secondary Structure**.

> **Hairpin screen versus folding.** The hairpin screen lists stretches
> that could pair with each other (stems of 4-7 nt, loops of 3-12 nt). It
> knows no energies, and most of its candidates are not stable on their
> own. Entry 6 computes the minimum-free-energy structure with ViennaRNA
> (Turner 2004 parameters), which is the better guide. Both describe one
> molecule in isolation, not the RNA in a cell.

### Step 4: Protein analysis

Choose **4**. The protein is translated from the selected ORF with its
genetic code, and the first codon is read as methionine (the initiator),
also when it is GTG or TTG. Entries 1-5 show the overview, the
translation, the amino-acid composition, physicochemical properties and
a complete report.

> **What the numbers mean.**
> - *Molecular weight*: the average mass of the chain.
> - *Isoelectric point (pI)*: the pH at which the protein carries no net
>   charge. It is computed from the sequence; modifications and folding
>   shift the real value.
> - *GRAVY*: the average Kyte-Doolittle hydropathy; positive values
>   suggest a hydrophobic, possibly membrane-associated protein.
> - *Instability index* (Guruprasad et al. 1990): above 40 a protein is
>   predicted to be unstable in a test tube. It is an empirical rule from
>   1990 data, not a measurement.
> - *Extinction coefficient*: absorbance at 280 nm, used to measure
>   protein concentration; it assumes the protein in water.

### Step 5: Restriction analysis and cloning

Choose **5**, then **1 Restriction Analysis**. Enter `EcoRI,HindIII` at
the enzyme prompt and keep the topology (circular). The digest cuts
pUC19 at the two ends of its multiple cloning site: two fragments of 51
and 2,635 bp. **3 Digest & Gel Analysis** also lists the predicted gel
bands, and Type IIS enzymes (BsaI, BsmBI, BbsI, SapI) for Golden Gate
assembly.

> **Sticky and blunt ends.** EcoRI cuts G^AATTC and leaves a 5' overhang
> AATT; two ends with matching overhangs can be ligated. Blunt cutters
> such as EcoRV leave no overhang. Type IIS enzymes cut outside their
> recognition site, so the overhang can be chosen freely.

**2 Cloning Analysis** combines the vector with an insert from the same
file. The automatic planner ranks strategies (single enzyme, two enzymes,
compatible ends) and explains each one. A strategy is a sequence-level
plan: ligation efficiency, methylation blocking and star activity are
not modeled.

### Step 6: Primer design

Choose **6**, then **1 Design Primers**, and give a target region, for
example 396 to 447 (the multiple cloning site). Synapse searches forward
and reverse primers around it and ranks pairs.

> **Melting temperature.** Synapse computes Tm with the nearest-neighbor
> model (SantaLucia 1998) at 50 mM Na⁺, 1.5 mM Mg²⁺, 0.2 mM dNTP and
> 250 nM primer; the conditions can be changed. Its values agree with
> Primer3 within 0.003 °C. A predicted Tm is not an annealing
> temperature: determine that with a gradient PCR.

Hairpin and dimer checks are screens (length of the longest
complementary stretch), not free energies. **2 Design Primers Around a
Restriction Site** designs pairs whose amplicon contains one chosen site
and predicts its digest. **3 Single Primer Analysis** checks a primer you
already have.

### Step 7: Report

Choose **P** to see what the session has stored (sequence, selected ORF,
RNA, protein, digest, primers) with the parameters of each step. **R**
exports everything as JSON (for programs) or Markdown (for people). The
report records the input's SHA-256 checksum, every parameter, the
genetic code and the Synapse version, so the analysis can be repeated.

## 6. Module reference

### Sequence Overview (1)

Length, composition, GC/AT content and skews, molecular weight, CpG
statistics and islands, homopolymers, ambiguous IUPAC symbols,
complement and reverse complement, GC and complexity profiles. Works for
DNA; RNA is analyzed in module 3.

### ORF Analysis (2)

Six-frame ORF search with any of the 27 NCBI genetic codes, start-codon
policies (ATG only, ATG/GTG/TTG, or the table's own starts), minimum
length, nested ORFs, partial ORFs at the ends of linear sequences and
ORFs across the origin of circular ones. Ranking is by length (default)
or by a heuristic score; the score is a convenience, not a probability.
*Compare ORFs* and *ORF Relationships* show overlaps and nested ORFs.

### RNA Analysis (3)

Transcription of the selected ORF or the whole sequence, composition,
codon usage and the three-frame translation in the ORF's genetic code,
Kozak context, hairpin screen, and the optional MFE structure.

### Protein Analysis (4)

From a selected ORF or a protein record: translation (ambiguous codons
give `?`), mass, pI, GRAVY, instability index, extinction coefficient,
composition, hydropathy profile, and transmembrane and signal-peptide
screens. The screens use the Kyte-Doolittle scale and are not
predictions from TMHMM or SignalP.

### Restriction & Cloning (5)

22 Type II enzymes and 4 Type IIS enzymes, IUPAC-aware site search (a
site with an N in the recognition sequence is reported as POSSIBLE, not
counted as a cut), linear and circular digests, fragment table, gel
bands, end compatibility, and the cloning planner. Enzyme data:
REBASE-checked (see METHODS.md).

### Primer Design (6)

Pair design with nearest-neighbor Tm, GC clamp, hairpin, self- and
cross-dimer screens, amplicon limits and specificity against the loaded
sequence. Three Tm models can be chosen (nearest-neighbor, Wallace,
salt-adjusted).

### Mutation & Variant Analysis (7)

**Mutation Analysis** aligns a reference and a variant sequence (for
example records 1 and 3 of the example file), lists every substitution,
insertion and deletion, and predicts coding consequences in the HGVS
notation (`c.7C>T`, `p.(R3C)`, `c.6dup`). The protein change is written
in parentheses because it is predicted from DNA. Position 1 must be the
first base of the coding sequence for c. numbering to equal HGVS c.
numbering.

> **Reading HGVS.** `c.215C>G` means: in the coding sequence, base 215,
> C replaced by G. `p.(Pro72Arg)` or `p.(P72R)`: proline 72 predicted to
> become arginine. `del`, `dup`, `ins` and `delins` describe deletions,
> duplications, insertions and replacements; in a repeat the change is
> written at its most 3' position (the HGVS 3' rule).

**Variant Analysis** reads VCF files, normalizes variants (left
alignment and trimming, as bcftools norm does), filters, computes
statistics such as the transition/transversion ratio, compares variant
sets and exports VCF or TSV.

### CRISPR Guide Design (8)

Guide candidates for SpCas9, SpCas9-NG, SaCas9 and Cas12a, with modeled
cut sites, base- and prime-editing candidates and an off-target search in
the sequences you supply (README.md). Limitations: no CFD score, no
genome-wide off-target search; the module is marked partially
implemented. Validate guides experimentally.

### BLAST Search (9)

Searches NCBI's servers or a local BLAST+ database.

- **Remote** searches send the query to NCBI; Synapse asks before it
  does. Do not submit unpublished or confidential sequences. NCBI can
  take several minutes; a search can be resumed later with its request
  ID (RID).
- **Local** searches need BLAST+ and a database. The database manager
  (entry 4) downloads NCBI databases and shows their size, release and
  your free disk space first; swissprot (about 0.4 GB) is a good start for
  proteins.
- Results show identity, query coverage, E-value and an evidence label.
  The **E-value** is the number of hits this good expected by chance in a
  database of this size: smaller is better, and 0 means "too small to
  represent". Similarity does not prove the same function.

### Literature Search (10)

Searches PubMed, Europe PMC, OpenAlex and Crossref together, removes
duplicates, ranks by a transparent relevance score and exports JSON,
CSV, text, RIS or BibTeX. Rankings are a starting point for reading, not a judgement
of quality.

### AI Experiment Planner (11)

Coming Soon. Planned, not available.

## 7. Reports, exports and reproducibility

Every analysis result carries:

- the input's SHA-256 checksum, so you can show which sequence was
  analyzed;
- every parameter, including defaults you did not change;
- the genetic code, the database and enzyme-table versions where they
  apply;
- warnings and limitations as text.

To repeat an analysis, load the same file (same checksum) and use the
parameters from the report. The JSON report can be read by other
programs; the Markdown report is for lab notebooks.

## 8. Using Synapse from Python

Every analysis is available without the menus:

```python
from core import api

record = api.load_records("data/examples/pUC19.fasta")[0]

orfs = api.find_orfs(record, genetic_code=11)
longest = orfs.data["orfs"][0]
print(longest["strand"], longest["start"] + 1, longest["stop"] + 1)   # - 1626 2486

digest = api.restriction_digest(record, ["EcoRI", "HindIII"])
print([fragment["length"] for fragment in digest.data["fragments"]])  # [51, 2635]

print(orfs.parameters)    # every parameter used
print(orfs.limitations)   # what the result does not show
```

Each function returns an `AnalysisResult` with `data`, `parameters`,
`warnings`, `limitations` and the input checksum. Returned coordinates
are 0-based (add 1 for the screen convention). The function list is in
`core/api.py`; the conventions are in METHODS.md.

## 9. How far the results can be trusted

Synapse distinguishes three levels (docs/CHANGES.md):

1. **Software correctness**: the code does what it specifies (about 700
   automated tests).
2. **Method conformity**: results agree with the published method or an
   independent tool.
3. **Biological validity**: confirmation by experiment. Synapse claims
   this for nothing.

What has been compared (details and dates in METHODS.md and
[BENCHMARKS.md](BENCHMARKS.md)):

| Area | Reference | Result |
|---|---|---|
| ORFs | NCBI ORFfinder (pUC19, 4 settings) | 69 of 69 identical |
| ORFs | RefSeq gene annotation of *E. coli* K-12 (4,294 genes) | every gene of at least 90 nt found at its stop codon |
| Genetic codes | NCBI gc.prt, 27 tables | identical |
| Primer Tm | Primer3 (2,916 comparisons) | within 0.003 °C |
| Restriction | REBASE, Biopython (Type IIS cut positions) | identical |
| Variant normalization | bcftools norm test vectors | 29 of 29 |
| HGVS | Mutalyzer; 1,500 ClinVar variants | c. 1,489 identical, 11 equivalent; p. all 1,186 |
| Protein | ExPASy ProtParam, Biopython | 399 of 400 instability weights; pI within 1e-4 |
| RNA structure | ViennaRNA RNAfold | MFE view uses ViennaRNA itself |
| BLAST | NCBI BLAST+ with refseq_rna; NCBI web service | expected hits, E-values and strands |

What has not been checked: experimental behavior of any primer, guide or
construct; genome-wide specificity; clinical interpretation.

## 10. Troubleshooting

| Problem | Cause and remedy |
|---|---|
| `ERROR: Invalid DNA sequence` | The record holds symbols that are not IUPAC DNA (digits, `-`, `*`). Clean the file or load it as protein. |
| Synapse asks whether a record is DNA or protein | The letters fit both. Answer the question; GenBank files avoid it. |
| "Secondary Structure ... not available" | ViennaRNA is not installed, or Windows Smart App Control blocks it. |
| A remote BLAST takes long | NCBI queues searches. Note the RID and resume later (BLAST menu, entry 2). |
| Local BLAST finds no database | Set `SYNAPSE_BLAST_DB_ROOT` to the folder that holds the databases, or download one in the database manager. |
| A whole genome is slow | A 4.6 Mb bacterial genome takes about 10 s for ORFs, about 25 s for a 22-enzyme digest and about 1 min for the overview (BENCHMARKS.md). |
| Strange characters in the terminal | Use Windows Terminal or PowerShell with UTF-8. |

Log files in `logs/` record what Synapse did; include them when you
report a problem (see CONTRIBUTING.md).

## 11. Glossary

| Term | Meaning |
|---|---|
| Amplicon | The DNA stretch a primer pair amplifies |
| Codon | Three bases that encode one amino acid or a stop |
| CDS | Coding sequence: the part of a gene translated into protein |
| E-value | Expected number of chance hits at least this good (BLAST) |
| Genetic code | The table assigning codons to amino acids; NCBI lists 27 |
| HGVS | Standard notation for sequence variants (den Dunnen et al. 2016) |
| IUPAC code | Letters for ambiguous bases: R = A/G, Y = C/T, N = any |
| MFE | Minimum free energy: the most stable predicted RNA structure |
| ORF | Open reading frame: start codon to in-frame stop codon |
| Overhang | Single-stranded end left by a staggered restriction cut |
| Reading frame | One of three ways to divide a strand into codons |
| Reverse complement | The opposite strand, read 5' to 3' |
| Tm | Melting temperature: half of the duplexes are dissociated |
| VCF | Variant Call Format: the standard file for sequence variants |

## References

- den Dunnen JT et al. (2016) HGVS recommendations for the description of
  sequence variants: 2016 update. *Hum Mutat* 37:564-569.
- Gardiner-Garden M, Frommer M (1987) CpG islands in vertebrate genomes.
  *J Mol Biol* 196:261-282.
- Guruprasad K, Reddy BVB, Pandit MW (1990) Correlation between stability
  of a protein and its dipeptide composition. *Protein Eng* 4:155-161.
- Kyte J, Doolittle RF (1982) A simple method for displaying the
  hydropathic character of a protein. *J Mol Biol* 157:105-132.
- Lorenz R et al. (2011) ViennaRNA Package 2.0. *Algorithms Mol Biol* 6:26.
- SantaLucia J Jr (1998) A unified view of polymer, dumbbell, and
  oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS*
  95:1460-1465.
