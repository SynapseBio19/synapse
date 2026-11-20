# Module splitting tools (Phase 7.1 and Phase 8)

These scripts moved `analyses/blast.py`, `primer.py`, `mutation.py` and
`cloning.py` into the internal packages `analyses/_<name>/`. They are kept
for the next modules (`docs/ARCHITECTURE.md` §8).

- `split_module.py <plan> [--check]` cuts the module named in the plan at
  statement boundaries.
  - Each top-level statement, with the comments before it, goes to the
    module its names are assigned to. A statement that defines no name
    goes with the module that owns the names it uses; `if __name__` blocks
    and the plan's `FACADE_NAMES` (such as `__all__`) stay in the public
    module.
  - The imports of every new module are generated from the names it uses.
    Several dotted imports under one name (`import urllib.error`,
    `import urllib.parse`) are all kept (fixed in Phase 10; before, only
    one was generated).
  - With `--check` it only prints the module graph and stops on a cycle.
  - It writes the package with LF line endings and rewrites the public
    module as docstring + re-exports + entry block, keeping that file's
    line endings.
- `verify_move.py <original> <public> <package_dir> [old=>new ...]`
  checks that every top-level statement of the original occurs verbatim
  in exactly one of the new files. The optional `old=>new` pairs are the
  plan's `REWRITES`.
- `plans/` holds the eleven plans as they were run (the last one,
  `literature.py`, in Phase 10). They are records: their `SOURCE` modules
  have since been split, so running them again would fail.
- A module whose tests replace helpers on it can still be split: record
  its characterization with the network or the operating system replaced
  below every Synapse function (the literature section replaces
  `urllib.request.urlopen`), so the snapshot reaches the code after the
  move as well.

**The procedure.** Each step was one commit, in this order:

1. Characterization snapshot (`tests/characterization.py`) and namespace
   list (`tests/test_module_namespaces.py`), recorded from the unsplit
   module.
2. `--check` until the graph has no cycle.
3. Split, `verify_move.py`, then the full suite on 3.14 and 3.10 and the
   golden master.
4. Search the tests for `patch.object(<module>, ...)`,
   `<module>.NAME = ...` and reads of module state. After a split such a
   test no longer reaches the code but usually still passes. Point it at
   the implementing module.
5. Add the package to `packages` in `pyproject.toml`, to the
   independence test and to the CI import check.
