## What changes and why

<!-- One logical change. For a scientific change: the reason, the reference and the evidence. -->

## Checklist (CONTRIBUTING.md)

- [ ] `python -m unittest discover -s tests -t .` passes (Python 3.10 and a current version)
- [ ] `python -m unittest -v tests.test_golden_master` passes, or `docs/CHANGES.md` explains each changed snapshot
- [ ] `python -m ruff check .` and `python scripts/ci/pyright_gate.py` pass
- [ ] A test fails on the old code and passes with the change (for a fix)
- [ ] `docs/METHODS.md`, `docs/FEATURE_STATUS.md` and `docs/CHANGES.md` describe the new behavior
- [ ] No fabricated reference values; every external value names its source, version and date
- [ ] Heuristic scores are labeled as heuristics
