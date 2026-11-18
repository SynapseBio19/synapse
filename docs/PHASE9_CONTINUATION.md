# Phase 9 continuation record (in progress)

**Phase 9: completion under the owner's mandate.** Started 2026-11-19.
The owner wants every defect fixed and the project completed. Only the
organizational steps (GitHub, licence) should remain for the owner.
CRISPR and the AI Experiment Planner stay excluded.

**Read first:**
- `docs/CHANGES.md`, section "Phase 9";
- `docs/ARCHITECTURE.md`;
- `docs/PHASE8_REPORT.md`.

## Environment (owner's machine)

| Item | State |
|---|---|
| Smart App Control | **off** since 2026-11-19, so primer3-py, ViennaRNA and `pyright.exe` load |
| BLAST databases | `C:\SynapseData\blast_db`, with `SYNAPSE_BLAST_DB_ROOT` set in the user environment |
| `refseq_rna` | VERIFIED: 83,968,556 sequences, dated 2026-09-30; NCBI has a release from 2026-10-30 |
| Python | 3.14.6 in `.venv`; 3.10.11 in `.tools/python-3.10.11` (`baseline.sh`-style runs) |

The test suite clears the `SYNAPSE_*` variables itself (`tests/__init__.py`).

## Done in Phase 9 (commits)

| Package | Commit | Content |
|---|---|---|
| A | `20f1dd8` | All 27 NCBI genetic codes, checked against `gc.prt` |
| B | `1fe380f` | Differential test of the ORF finder and RNA translation in all codes |
| C | `ec459e2` | Size, release and free-space check before an NCBI database download |
| D | `3ee20a6` | `synapse` command; `core/paths.py`; no user files in site-packages; clean exit on end of input; three latent defects |

## Next actions, in this order

1. **Split the remaining large modules** with `scripts/refactor/`
   (procedure in its README), one commit each: `protein`,
   `variant_analysis`, `sequence`, `orf`, `restriction`, `rna`,
   `literature`.
   - Before each split: a characterization section (`tests/characterization.py`)
     and the namespace list (`tests/test_module_namespaces.py`).
   - After each split: search the tests for patches and module-state reads
     of the public module.
2. Consolidate the duplicated sequence helpers on `core.seqtools`.
3. Final validation:
   - suite on 3.14 and 3.10, also with `SYNAPSE_BLAST_DB_ROOT` set;
   - golden master; Ruff; pyright gate;
   - the CI import step run exactly as written against a fresh install;
   - optionally `SYNAPSE_REFSEQ_RNA_TESTS=1` (17 minutes).
4. The Phase 9 report.
