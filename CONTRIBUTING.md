# Contributing to Synapse

Synapse is research and teaching software whose results must be
traceable. These rules come from the project's history
(`docs/PROJECT_HISTORY_AND_DECISIONS.md`) and are checked by the test
suite where possible.

**Licence.** Synapse is licensed under the GNU Affero General Public
License, version 3 or later (`LICENSE`). Contributions are accepted under
the same licence. The copyright holder also offers Synapse under other
terms; before an outside contribution is merged, the contributor may be
asked to agree that it can be included under those terms as well
(`docs/LICENSING_AND_GOVERNANCE.md`, section 6).

## Setup

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -e ".[dev,reference,structure]"
```

Use `.venv/bin/python` on Linux and macOS. The application itself needs
only the Python standard library (3.10 or newer); the extras add the
reference packages used by the tests and the development tools. NCBI
BLAST+ is optional.

## Before every commit

```bash
python -m unittest discover -s tests -t .
python -m unittest -v tests.test_golden_master
python -m ruff check .
python scripts/ci/pyright_gate.py
```

The CI (`.github/workflows/ci.yml`) runs the same commands on Python
3.10-3.14.

## Rules

- **Science before cosmetics.** Do not change a scientific result without
  a reason, a reference and a test. A test that fails on the old code is
  the expected evidence for a fix.
- **Golden master.** Snapshots in `tests/golden/` change only together
  with an entry in `docs/CHANGES.md` that explains the scientific reason
  (`SYNAPSE_UPDATE_GOLDEN=1`).
- **No fabricated data.** Reference values, citations and validation
  results must come from a source that is named, with its version and
  date (`docs/TOOL_AND_DATA_LEDGER.md`).
- **Heuristics are labeled.** A score that is not a validated probability
  is called a heuristic in code, output and documentation.
- **Layers.** Science lives in `analyses/` and `core/` and never imports
  `ui`; dialogs live in `ui/` and draw through `ui/style.py`
  (`docs/TERMINAL_UI_STANDARD.md`). Messages use the ERROR, WARNING, NOTE
  and OK prefixes.
- **Language.** Code, comments, output and documentation are in English,
  in American spelling (checked by the suite).
- **Optional dependencies stay optional.** `dependencies = []` in
  `pyproject.toml`; reference packages are imported inside the tests or
  behind a check that reports their absence.
- **Documentation follows the code.** Update `docs/FEATURE_STATUS.md`,
  `docs/METHODS.md` and `docs/CHANGES.md` in the same commit. A feature
  is "Execution-tested" only if a test runs it and checks its output.
- **History.** One logical change per commit, with the evidence in the
  message. No history rewrites or force-pushes.
- **Owner decisions** (binding until the owner changes them):
  - The AI Experiment Planner stays "Coming Soon".
  - CRISPR-related work is paused.
