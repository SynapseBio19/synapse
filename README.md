# Synapse

Sequence-centric DNA analysis workbench (command-line, Python standard library only).

Status: research and teaching software under active correction. Results are
computational predictions; none is experimentally validated, and Synapse is
not a substitute for established, validated tools. See `docs/METHODS.md` for
methods, defaults, references and validation status, and `docs/CHANGES.md`
for every change to scientific behavior. `docs/AUDIT_STATUS.md` records what
the 2026-10 audit asked for and what has been implemented, tested, deferred
or rejected.

Requirements: Python 3.10 or newer; no third-party packages.

## Run

```
python main.py
```

Load a FASTA or GenBank file (Enter loads `examples/dna.fasta`), choose a
sequence and an analysis. `[P]` shows stored results, `[R]` exports a
reproducible report (JSON or Markdown).

## Use without the menu

```python
from core import api

records = api.load_records("examples/dna.fasta")
orfs = api.find_orfs(records[0])                 # AnalysisResult
print(orfs.parameters, orfs.data["candidate_count"])
digest = api.restriction_digest(records[0], ["EcoRI", "BamHI"])
```

Every function returns an `AnalysisResult` with parameters, input SHA-256,
warnings and limitations. No function asks for keyboard input.

## Layout

- `analyses/` – scientific modules (each also has its terminal menu)
- `core/` – shared models, genetic codes, pipeline state, report, API
- `fileio/` – FASTA/GenBank loading
- `data/` – enzyme and genetic-code data
- `tests/` – unit, regression (golden master) and end-to-end tests

## Test

```
python -m unittest discover -s tests -t .
```

The suite blocks all network access. Golden-master snapshots in
`tests/golden/` change only with a documented reason
(`SYNAPSE_UPDATE_GOLDEN=1`). BLAST tests that need BLAST+ or local
databases (`data/blast_db/`) are not part of the default run.
