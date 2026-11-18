# Legacy BLAST record – test evidence, do not use as a result

Moved here from `blast_results/` on 2026-11-18 (Phase 8): the file is
test data, not a result directory.

`TEST_8_REFERENCE_blast.json` is the stored output of one remote NCBI
BLAST run (RID `CD3MMW95016`, 2026-10-07, blastn against `refseq_rna`),
written by Synapse **before** the P0-1 correction. It is kept unchanged
as the record of that run and of the defect, not as a valid result
(audit D7, `docs/PHASE6_INITIAL_AUDIT.md`):

- `best_coverage_percent` is 100.0 for a 32 bp HSP on a 99 bp query. The
  corrected definition (union of HSP query intervals / query length)
  gives 32.32 %, which `tests/test_p0_regressions.py` checks with these
  values.
- The hit has E = 7.45 and no evidence band; current Synapse assigns
  bands only at E ≤ 1e‑3 and reports this hit as `NOT_SIGNIFICANT`.
- The HSP keys (`query_from`/`query_to`) are an older schema; current
  code (`query_start`/`query_end`) cannot re-read the file.

`tests/test_p0_regressions.py` (`test_stored_test8_result_recomputed`)
reads it to show the corrected coverage; no application code does. (This
note said "no test reads this file" until 2026-11-18, which was wrong.) A
current run is reproduced by the BLAST tests in `tests/test_phase6_tier2.py`
and `tests/test_phase6_sequence_types.py` (local BLAST+).
