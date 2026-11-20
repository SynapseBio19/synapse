# Synapse

A sequence-centric DNA, RNA and protein analysis workbench for the
terminal, written in Python with the standard library only.

## Purpose and status

Synapse brings the routine sequence analyses of a molecular biology
project into one session: load a sequence once, analyze it from several
angles, keep every result with its parameters, and export a reproducible
report. It is designed to be transparent about what each analysis
computes and what it does not.

New users start with the user guide, `docs/USER_GUIDE.md`: concepts,
every menu, and a tutorial that takes the plasmid pUC19 from loading to
report.

Synapse is research and teaching software under active correction.
**Its results are computational; none has been validated
experimentally**, and it is not a substitute for established, validated
tools. Where a module has been compared with an independent reference,
`docs/METHODS.md` says which one and with what result.

## Analysis modules

| Menu | Module | What it does |
|---|---|---|
| [1] | Sequence overview | Composition, GC/AT content and skew, ambiguity codes, homopolymers, CpG, molecular mass |
| [2] | ORF analysis | ORF search on both strands with a selectable NCBI genetic code, start-codon policy and minimum length; ORFs across the origin of circular sequences and optional 5'-partial ORFs; candidate comparison and selection |
| [3] | RNA analysis | Transcription of the selected ORF or sequence, codon view with the ORF's genetic code, hairpin screen (Watson-Crick, G·U counted separately), sequence motifs; optional MFE secondary structure with ViennaRNA |
| [4] | Protein analysis | Translation and sequence-derived protein descriptors; protein records open this module directly |
| [5] | Restriction and cloning | Restriction map and digest (linear or circular) with 22 palindromic and 4 Type IIS enzymes, enzyme pairs, restriction-cloning strategies and an automatic planner |
| [6] | Primer design | Primer pairs with nearest-neighbor Tm and salt correction, structural screens, ranking; primers around a chosen restriction site with the predicted digest of each amplicon |
| [7] | Mutation and variant analysis | Sequence comparison with coding consequences in HGVS-style notation; VCF parsing, reference checks, normalization, filtering |
| [8] | CRISPR guide design | Guides for SpCas9, SpCas9-NG, SaCas9 and Cas12a, modeled cut sites, base and prime editing candidates, off-target search in supplied sequences |
| [9] | BLAST search | Local BLAST+ or NCBI BLAST, with resumable request IDs and an evidence summary per hit |
| [10] | Literature search | PubMed, Europe PMC, OpenAlex and Crossref, deduplicated |
| [11] | AI Experiment Planner | **Coming Soon.** Planned, not yet available; the entry runs nothing. |

`[P]` lists the stored results, `[R]` exports the session report (JSON or
Markdown), `[S]` selects another sequence. `docs/FEATURE_STATUS.md` lists
every function with its implementation, test and validation status.

## Installation

Synapse needs Python 3.10 or newer and nothing else to run.

```bash
git clone <repository> synapse
cd synapse
python main.py
```

Or install it as a command (the package is called `synapse-gen`, the
command `synapse`):

```bash
python -m pip install .
synapse
```

An installed Synapse keeps logs, exports and local BLAST databases in a
per-user folder (`%LOCALAPPDATA%\Synapse` on Windows,
`~/.local/share/synapse` elsewhere); run from the source folder it uses
the project folder, as before (see Configuration).

Optional, for the reference comparisons in the test suite:

```bash
python -m pip install -e ".[reference]"
```

Optional, for the energy-based RNA secondary structure (RNA menu entry 6):

```bash
python -m pip install -e ".[structure]"
```

Development setup in a project-local virtual environment (no system-wide
installation; the commands are those the CI runs):

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -e ".[dev,reference,structure]"
```

On Linux and macOS the interpreter is `.venv/bin/python`. NCBI BLAST+ is a
separate download from NCBI and is needed only for local BLAST and the
tests that use it.

## Dependencies

| Kind | What | Needed for |
|---|---|---|
| Runtime | Python 3.10+ standard library | everything |
| Optional tool | NCBI BLAST+ (`blastn`, `blastp`, `makeblastdb`, `blastdbcmd`) | local BLAST searches |
| Optional package | ViennaRNA (`[structure]` extra) | MFE secondary structure in the RNA module; without it the view says it is not available |
| Optional network | NCBI BLAST and the four literature services | remote BLAST and literature search, only when the user starts one |
| Test reference | Biopython, primer3-py (`[reference]` extra) | independent cross-checks in the tests |
| Development | pyright, Ruff, coverage.py (`[dev]` extra) | type checking against Python 3.10, linting, line coverage |

The scientific core uses the standard library only, so that every
calculation is visible in this repository.

## Usage

### Terminal

```bash
python main.py
```

Load a FASTA or GenBank file (Enter loads the example `data/examples/dna.fasta`;
`data/examples/pUC19.fasta` is the circular cloning vector pUC19), choose a
record, then an analysis. A record's type (DNA, RNA or protein) is taken
from the file when it states one and otherwise from its letters; when the
letters cannot decide, Synapse asks rather than guesses. Positions shown
in the menus are 1-based and inclusive.

### Python

Every analysis is also available without the menus, through `core.api`:

```python
from core import api

records = api.load_records("data/examples/dna.fasta")
orfs = api.find_orfs(records[0])
print(orfs.parameters, orfs.data["candidate_count"])
digest = api.restriction_digest(records[0], ["EcoRI", "BamHI"])
```

Each function returns an `AnalysisResult` holding the data, the
parameters, the input's SHA-256, warnings and limitations. No API
function reads from the keyboard. Coordinates in returned data follow
the conventions documented per structure in `docs/METHODS.md`.

## Configuration

Analysis parameters are set per run, in the module menus or as API
arguments; their defaults and sources are listed in `docs/METHODS.md`.
The environment variables below configure the session:

| Variable | Effect |
|---|---|
| `SYNAPSE_HOME` | Synapse's folder for logs, exports and BLAST databases (default: the project folder from source, a per-user folder when installed) |
| `SYNAPSE_LOG_DIR` | Directory for log files (default `logs/` in Synapse's folder) |
| `SYNAPSE_EXPORT_DIR` | Directory for CRISPR exports (default `exports/crispr/`; a configured directory gets a `crispr/` subdirectory) |
| `SYNAPSE_BLAST_BIN`, `NCBI_BLAST_BIN`, `BLAST_HOME`, `NCBI_BLAST_HOME` | Where to look for BLAST+ executables |
| `SYNAPSE_BLAST_DB_ROOT` | Root directory of local BLAST databases (default `data/blast_db` from source, `blast_db` in the per-user folder when installed). NCBI databases are large (swissprot 0.4 GB, refseq_rna about 81 GB, core_nt about 308 GB); keep them outside cloud-synchronized folders such as OneDrive. The database manager shows the size and the free space before a download |
| `SYNAPSE_NCBI_EMAIL`, `SYNAPSE_NCBI_TOOL`, `SYNAPSE_NCBI_API_KEY` | Identification sent to NCBI, as NCBI asks of API users |
| `SYNAPSE_OPENALEX_EMAIL`, `SYNAPSE_CROSSREF_EMAIL` | Contact address for OpenAlex and Crossref |

## Architecture

```
main.py      session: file loading, record selection, main menu
analyses/    one module per analysis: the science, without any terminal code;
             the larger modules keep their code in internal packages
             analyses/_<name>/ (docs/ARCHITECTURE.md)
core/        data model, sequence types, pipeline state, report, API, text layout
fileio/      FASTA and GenBank reading
ui/          terminal interface: shared presentation (ui/style.py), main menu,
             report screens and one dialog module per analysis (ui/<module>_menu.py)
data/        enzyme definitions, genetic codes, amino-acid data
tests/       unit, golden-master, integration and end-to-end tests
docs/        methods, changes, audit and feature status, interface and style standards
```

No analysis module imports `ui/`: every analysis runs through `core.api`
without the terminal interface, which a test checks with `ui` made
unimportable.

`core/pipeline.py` keeps the session's sequence, the selected ORF and
every stored result, and invalidates RNA and protein results when the
selected ORF changes. `core/report.py` turns the stored results into the
exported report, with software version, input provenance and each
result's parameters, warnings and limitations.

`LEGACY/` (pre-project scripts) and `Synapse_Website/` (a static
placeholder page) are not part of the software and are excluded from
packaging and type checking.

## Testing

```bash
python -m unittest discover -s tests -t .
```

The suite has run on Python 3.14 and on Python 3.10.11, the oldest
supported version (`docs/TOOL_AND_DATA_LEDGER.md`). Line coverage is
recorded in `docs/TEST_COVERAGE.md`.

Checks the CI runs (`.github/workflows/ci.yml`):

```bash
python -m unittest -v tests.test_golden_master
python -m ruff check .
python scripts/ci/pyright_gate.py
```

`scripts/ci/pyright_gate.py` fails on any pyright error outside a short,
documented allow-list (the findings in the CRISPR module, which was not
changed in Phase 7). `scripts/ci/check_skips.py` lists skipped tests with
their reasons. The workflow tests Python 3.10 to 3.14 on Linux and 3.14 on
Windows, plus a job with every optional package. It has not yet run on a
hosted service: the repository has no remote.

The suite blocks all network access, so no test can depend on an external
service. Golden-master snapshots in `tests/golden/` protect every
module's output and change only with a documented reason
(`SYNAPSE_UPDATE_GOLDEN=1`, recorded in `docs/CHANGES.md`). Tests that
need BLAST+ or the reference libraries are skipped, with the reason
shown, when those are not installed.

## Limitations

The full list, per module, is in `docs/METHODS.md`. The most important:

- No result is experimentally validated. Thermodynamic values, guide
  scores, signal-peptide and structure screens are model predictions or
  heuristics, and are labeled as such.
- Off-target search covers only the sequences supplied, never a genome.
- ORF calls were compared with NCBI ORFfinder on pUC19 (69 ORFs
  identical; Synapse counts the stop codon in the minimum length,
  ORFfinder does not) and with the RefSeq annotation of *E. coli* K-12
  (every annotated gene of at least 90 nt found at its stop codon). ORFs
  across the origin were not part of either comparison. An ORF search is
  not a gene predictor.
- The RNA hairpin screen has no energy model of its own; the optional
  MFE view (ViennaRNA) supplies one.
- HGVS descriptions number the supplied sequence; they equal HGVS c./p.
  numbering only if that sequence is the CDS. On 1,500 ClinVar variants,
  11 long deletion-insertions are written as an equivalent multi-change
  allele instead of ClinVar's single delins (docs/BENCHMARKS.md).
- Restriction analysis covers 22 palindromic and 4 Type IIS enzymes;
  methylation blocking and star activity are not modeled. Golden Gate
  assembly is not planned.
- The global alignment of the mutation module keeps one byte per matrix
  cell; similar sequences are aligned in an exact band (10 kb in well
  under a second), unrelated long sequences still take time proportional
  to the product of the lengths.
- Remote BLAST and literature search depend on external services, their
  availability and their rate limits.

## Licence

No licence has been chosen yet; until the owner decides, all rights are
reserved. The options and the decisions they need are in
`docs/LICENSING_AND_GOVERNANCE.md`. How to contribute: `CONTRIBUTING.md`.

## Documentation

- `docs/USER_GUIDE.md`: the user guide (concepts, menus, tutorial,
  interpretation, troubleshooting)
- `docs/METHODS.md`: algorithms, equations, defaults, references,
  assumptions, limitations and validation evidence
- `docs/BENCHMARKS.md`: genome-scale performance and the benchmarks on
  curated data (RefSeq gene models, ClinVar)
- `docs/RELEASE_NOTES.md`: the prepared release 0.3.0 and the steps before
  it
- `CITATION.cff`: how to cite Synapse
- `docs/FEATURE_STATUS.md`: every feature with its status
- `docs/CHANGES.md`: every change to behavior, with reasons
- `docs/AUDIT_STATUS.md`: the status of every audit recommendation
- `docs/STYLE_GUIDE.md` and `docs/TERMINAL_UI_STANDARD.md`: communication
  and interface standards
- `docs/TEST_COVERAGE.md` and `docs/TOOL_AND_DATA_LEDGER.md`: test
  coverage, and every external tool and dataset with source, version and
  checksum
- `docs/PHASE6_FINAL_REPORT.md`: what Phase 6 changed, what was verified
  and what remains
