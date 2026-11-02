# Synapse

A sequence-centric DNA, RNA and protein analysis workbench for the
terminal, written in Python with the standard library only.

## Purpose and status

Synapse brings the routine sequence analyses of a molecular biology
project into one session: load a sequence once, analyze it from several
angles, keep every result with its parameters, and export a reproducible
report. It is designed to be transparent about what each analysis
computes and what it does not.

Synapse is research and teaching software under active correction.
**Its results are computational; none has been validated
experimentally**, and it is not a substitute for established, validated
tools. Where a module has been compared with an independent reference,
`docs/METHODS.md` says which one and with what result.

## Analysis modules

| Menu | Module | What it does |
|---|---|---|
| [1] | Sequence overview | Composition, GC/AT content and skew, ambiguity codes, homopolymers, CpG, molecular mass |
| [2] | ORF analysis | ORF search on both strands with a selectable NCBI genetic code, start-codon policy and minimum length; candidate comparison and selection |
| [3] | RNA analysis | Transcription of the selected ORF or sequence, codon view, sequence motifs |
| [4] | Protein analysis | Translation and sequence-derived protein descriptors; protein records open this module directly |
| [5] | Restriction and cloning | Restriction map and digest (linear or circular), enzyme pairs, restriction-cloning strategies |
| [6] | Primer design | Primer pairs with nearest-neighbor Tm and salt correction, structural screens, ranking |
| [7] | Mutation and variant analysis | Sequence comparison with coding consequences; VCF parsing, reference checks, normalization, filtering |
| [8] | CRISPR guide design | Guides for SpCas9, SpCas9-NG, SaCas9 and Cas12a, modeled cut sites, base and prime editing candidates, off-target search in supplied sequences |
| [9] | BLAST search | Local BLAST+ or NCBI BLAST, with resumable request IDs and an evidence summary per hit |
| [10] | Literature search | PubMed, Europe PMC, OpenAlex and Crossref, deduplicated |
| [11] | Experiment planner | **Not implemented.** The entry says so. |

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

Optional, for the reference comparisons in the test suite:

```bash
python -m pip install -e ".[reference]"
```

## Dependencies

| Kind | What | Needed for |
|---|---|---|
| Runtime | Python 3.10+ standard library | everything |
| Optional tool | NCBI BLAST+ (`blastn`, `blastp`, `makeblastdb`, `blastdbcmd`) | local BLAST searches |
| Optional network | NCBI BLAST and the four literature services | remote BLAST and literature search, only when the user starts one |
| Test reference | Biopython, primer3-py (`[reference]` extra) | independent cross-checks in the tests |
| Development | pyright (`[dev]` extra) | type checking against Python 3.10 |

The scientific core uses the standard library only, so that every
calculation is visible in this repository.

## Usage

### Terminal

```bash
python main.py
```

Load a FASTA or GenBank file (Enter loads `examples/dna.fasta`), choose a
record, then an analysis. A record's type (DNA, RNA or protein) is taken
from the file when it states one and otherwise from its letters; when the
letters cannot decide, Synapse asks rather than guesses. Positions shown
in the menus are 1-based and inclusive.

### Python

Every analysis is also available without the menus, through `core.api`:

```python
from core import api

records = api.load_records("examples/dna.fasta")
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
| `SYNAPSE_LOG_DIR` | Directory for log files (default `logs/`) |
| `SYNAPSE_EXPORT_DIR` | Directory for CRISPR exports (default `exports/crispr/`; a configured directory gets a `crispr/` subdirectory) |
| `SYNAPSE_BLAST_BIN`, `NCBI_BLAST_BIN`, `BLAST_HOME`, `NCBI_BLAST_HOME` | Where to look for BLAST+ executables |
| `SYNAPSE_BLAST_DB_ROOT` | Root directory of local BLAST databases |
| `SYNAPSE_NCBI_EMAIL`, `SYNAPSE_NCBI_TOOL`, `SYNAPSE_NCBI_API_KEY` | Identification sent to NCBI, as NCBI asks of API users |
| `SYNAPSE_OPENALEX_EMAIL`, `SYNAPSE_CROSSREF_EMAIL` | Contact address for OpenAlex and Crossref |

## Architecture

```
main.py      session: file loading, record selection, main menu
analyses/    one module per analysis, each with its terminal menu
core/        data model, sequence types, pipeline state, report, API
fileio/      FASTA and GenBank reading
ui/          shared terminal presentation (ui/style.py), main menu, report screens
data/        enzyme definitions, genetic codes, amino-acid data, BLAST test database
tests/       unit, golden-master, integration and end-to-end tests
docs/        methods, changes, audit and feature status, interface and style standards
```

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
- ORF search does not report 5'-partial ORFs or ORFs spanning the origin
  of a circular sequence.
- Restriction analysis searches the forward strand, which is complete for
  the palindromic enzymes in the database; it holds no Type IIS enzyme.
  Methylation blocking and star activity are not modeled.
- The global sequence alignment in the mutation module needs memory in
  proportion to the product of the two lengths (about 4.5 GB for two
  10 kb sequences).
- Remote BLAST and literature search depend on external services, their
  availability and their rate limits.

## Documentation

- `docs/METHODS.md`: algorithms, equations, defaults, references,
  assumptions, limitations and validation evidence
- `docs/FEATURE_STATUS.md`: every feature with its status
- `docs/CHANGES.md`: every change to behavior, with reasons
- `docs/AUDIT_STATUS.md`: the status of every audit recommendation
- `docs/STYLE_GUIDE.md` and `docs/TERMINAL_UI_STANDARD.md`: communication
  and interface standards
- `docs/PHASE6_FINAL_REPORT.md`: what Phase 6 changed, what was verified
  and what remains
