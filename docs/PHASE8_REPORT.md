# Phase 8 report: consolidation under the owner's mandate

Date: 2026-11-18. Start: `3915d94` (end of Phase 7.1).

**Mandate** (owner, chat of 2026-11-18). In the copy, Claude may:
- restructure code;
- correct defects that are backed by evidence;
- delete dead code;
- rename internal names;
- drop the ReadOnly protection.

CRISPR and the AI Experiment Planner stay excluded for now. The name
**Synapse** stays, and the owner's domain is `synapse-gen.com`.

## 1. Summary

The four largest scientific modules now have the same structure: a thin
public module and an internal package with one module per responsibility.
- Phase 7.1: BLAST and primer.
- This phase: mutation and cloning.

Each module was frozen first (characterization snapshot and namespace
test) and then moved verbatim. Dead code is gone, the historical `_vN_`
names are replaced, and the package has a distribution name that can be
published.

| Metric | Phase 7.1 end | Phase 8 end |
|---|---:|---:|
| Tests | 601 | 608 |
| Failures | 0 | 0 |
| Coverage | 82.3 % | 83.7 % |
| Largest module in `analyses/` | `primer.py` was 5,182 lines; largest left: `mutation.py`, 4,670 | `protein.py`, 2,544 lines |
| Characterization scenarios | 78 | 182 |
| Pyright | 11 (crispr.py only) | 11 (crispr.py only) |
| Ruff | clean | clean |
| Golden master | 12/12 | 12/12 |

## 2. Packages

| Commit | Package | Content |
|---|---|---|
| `80ae423` | A | ReadOnly dropped (42 files). Deleted `print_protein_report` (599 lines), `print_primer_report` (39) and the root file `blast_test.json`. The legacy BLAST record moved to `tests/fixtures/legacy_blast/`, because a test reads it (its README had said none did). Two names that do not exist removed from `PRIMER_V3_PUBLIC_API`; a test now checks the list |
| `a7249fa` | B | 24 `_vN_` primer names made descriptive; public functions defined directly instead of bound by assignment |
| `6d4c216` | C | `mutation.py` (4,670 lines) -> `analyses/_mutation/` (11 modules), after a 79-scenario snapshot; float tolerance for the 3.10/3.14 `math.log2` difference |
| `4395a4a`, `962f7ee` | D | `cloning.py` (3,972 lines) -> `analyses/_cloning/` (7 modules), after a 25-scenario snapshot |
| `b3ba026` | E | pip distribution `synapse-gen`, homepage `https://synapse-gen.com`; the CI import check covers all four packages |
| final | F | `docs/ARCHITECTURE.md` (living document), `scripts/refactor/` (the splitting tools and plans), status documents, this report |

## 3. Defects and risks found

1. **Tests that would have passed without testing anything.** After a
   split, a test that patches or reads the public module no longer
   reaches the implementation, but it usually still passes. Three tests
   were affected:
   - **`test_alignment_memory`:** would have compared the banded alignment
     with itself.
   - **`test_planner_cache`:** would have counted no calls and still
     passed, and would have read a cache variable that never changes.

   All three now patch the implementing module. The counting test also
   asserts that it counted something.
2. **A self-test that checked `globals()`.** `run_cloning_self_test`
   verified the API through the globals of its own module. That holds only
   while the self-test lives in the public module. It now asks
   `analyses.cloning`.
3. **Wrong documentation.** `blast_results/README.md` claimed that no test
   read the legacy file. `PRIMER_V3_PUBLIC_API` listed two functions that
   `analyses.primer` does not have. Both corrected, the second with a
   test.
4. **Interpreter difference.** CPython 3.10.11 and 3.14.6 on Windows differ
   in the last digit of `math.log2`. The unchanged original code gives
   98.82341858018908 versus 98.82341858018907. Not a Synapse defect. The
   snapshot test tolerates a relative 1e-12, and only for floats.
5. **Environment.** Windows application control blocked ViennaRNA's DLL on
   2026-11-18 (primer3-py loaded the same day). Six tests skip with that
   reason; nothing was bypassed.

## 4. Verification (final state)

- **Suite:** 608 tests, 0 failures.
  - CPython 3.14.6: 7 skipped (benchmark, ViennaRNA).
  - CPython 3.10.11: 23 skipped (optional packages not installed there).
- **Snapshot:** all 182 characterization scenarios unchanged against
  their recordings.
- **Golden master:** 12/12.
- **Static checks:** Ruff clean; pyright 11 errors, all in the excluded
  `analyses/crispr.py`.
- **Clean non-editable install:** `pip show synapse-gen` reports the
  homepage. All four packages import from `site-packages`, and the BLAST
  and cloning self-tests pass there.
- **CI:** still never run on GitHub (no remote).

## 5. Not done, and why

| Item | Reason / next step |
|---|---|
| Licence | Legal choice of the owner; the options are in `docs/LICENSING_AND_GOVERNANCE.md`. "Ich überlasse alles dir" was not taken as choosing a licence on the owner's behalf |
| GitHub remote, running CI | Needs the owner's account |
| Installable `synapse` command | `main.py` reads `examples/` next to itself; needs packaged example data and `main` in a package |
| Splitting `protein`, `variant_analysis`, `sequence`, `orf`, `restriction`, `rna`, `literature` | Same procedure (`scripts/refactor/README.md`); each is one session-sized step |
| Duplicated helpers (`reverse_complement` 7x) | Needs per-function characterization of error texts first |
| Long dialog functions | Byte-identical menu tests exist; a mechanical, low-value step |
| CRISPR, AI Experiment Planner | Excluded by the owner |
| Scientific benchmarks on curated data, lab validation | Need data sets, downloads and, for the lab, people |
