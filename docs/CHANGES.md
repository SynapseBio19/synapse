# Synapse change log (scientific and behavioral changes)

Every change that alters a scientific result or a default is recorded here
with the reason, the old and new behavior, and the evidence. Golden-master
snapshots (`tests/golden/`) are regenerated only together with an entry here.

Verification levels used below:

- **Software correctness** – the code does what its specification says (unit tests).
- **Method conformity** – the result agrees with the published method or an
  independent reference implementation (known-answer / cross-validation tests).
- **Biological validity** – the result is confirmed by experiment or curated
  data. Synapse does not claim this for any output.

## Stage 0 – baseline (2026-10-29)

- Git repository initialised; tag `baseline-audit` = project exactly as audited.
- Added `tests/` (standard-library `unittest`, no extra dependencies):
  - `test_builtin_selftests.py` runs the 13 built-in module self-tests.
  - `test_golden_master.py` + `golden/*.json`: snapshots of current behavior
    of all 12 analysis modules.
  - `test_p0_regressions.py`: known-answer tests for the nine P0 audit
    findings, initially marked `expectedFailure` (20 tests reproduce the defects).
- Added `pyproject.toml` (Python >= 3.10, no runtime dependencies).
- No source file changed.

### Compatibility fix: `analyses/protein.py` (Python 3.10/3.11)

- Problem: seven f-strings in report-printing functions spanned several lines
  inside `{...}`. That syntax is only valid from Python 3.12 (PEP 701); on
  3.10/3.11 the whole protein module failed to import. Found by running the
  new test suite on Python 3.10.
- Change: the expressions were joined onto one line. Pure syntax change, no
  behavior change (golden master unchanged).
