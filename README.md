# Synapse

Sequence-centric DNA analysis workbench (command-line, Python standard library only).

Status: research and teaching software under active correction. Results are
computational predictions; none is experimentally validated. See
`docs/CHANGES.md` for every change to scientific behavior.

## Run

```
python main.py
```

## Test

```
python -m unittest discover -s tests -t .
```

The suite blocks all network access. BLAST tests that need BLAST+ or local
databases (`data/blast_db/`) are not part of the default run.
