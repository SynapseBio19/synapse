# Tool and data ledger

Every external tool, package and dataset that Phase 6 used for testing or
scientific validation, with its source, version, date, licence and
checksum (Master Prompt 1.3 "Traceability and reproducibility", 1.3.1,
11.8). Recorded on 2026-11-12. "Not recorded" means the information was
not captured when the resource was obtained and cannot be reconstructed
from the repository.

None of these is a runtime dependency: `pyproject.toml` declares
`dependencies = []`. The optional extras are `reference` (Biopython,
primer3-py) and `dev` (pytest, pyright).

## Interpreters

| Resource | Version | Source | Obtained | Licence | Checksum | Purpose |
|---|---|---|---|---|---|---|
| CPython (project `.venv`) | 3.14.6 (tags/v3.14.6:c63aec6, MSC v.1944, 64-bit) | python.org installer, per-user install | Before Phase 6; date not recorded | PSF-2.0 | not recorded | Development, the whole suite, pyright, coverage |
| CPython, Windows embeddable package | 3.10.11 (tags/v3.10.11:7d4cc5a, 64-bit) | <https://www.python.org/ftp/python/3.10.11/python-3.10.11-embed-amd64.zip> | 2026-11-12 | PSF-2.0 | MD5 `f1c0538b060e03cbb697ab3581cb73bc` (equal to the value on the python.org release page); SHA-256 `608619f8619075629c9c69f361352a0da6ed7e62f83a0e19c63e0ea32eb7629d` | Runtime test of the minimum supported version (`requires-python >= 3.10`) |

The 3.10 interpreter is unpacked in `.tools/python-3.10.11/` (git-ignored);
no system-wide installation was made. Its `python310._pth` was given the
project root as a third path entry so that the suite can import the
project. Result on 2026-11-12: `python -m unittest discover -s tests -t .`
ran 452 tests, all passing; 17 were skipped (the Biopython and primer3-py
comparisons, because those packages are not installed for 3.10; one test
that needs the 3.12 tokenizer; and the opt-in benchmark). The same suite
under 3.14.6 skips only the benchmark.

## Python packages (project `.venv`, from PyPI)

| Package | Version | Licence | Obtained | Used by |
|---|---|---|---|---|
| biopython | 1.88 | Biopython License Agreement (BSD-style) | Before Phase 6 (in `.venv` on 2026-10-30); exact date not recorded | `tests/test_external_reference.py`, `tests/test_phase6_tier1.py` (ORF translation), audit §3 |
| numpy | 2.5.4 | BSD-3-Clause and bundled licences | Dependency of biopython | – |
| primer3-py | 2.4.0 | GPL-2.0-or-later | 2026-11-02 | `tests/test_primer_thermodynamics.py`, `tests/test_primer_structure_reference.py` (used only as a test-time reference; Synapse neither imports nor distributes it at run time) |
| pyright | 1.1.414 | MIT | 2026-10-30 | Static analysis (`python -m pyright --pythonversion 3.10`) |
| nodejs-wheel-binaries | 24.19.0 | MIT | Dependency of pyright | – |
| coverage | 7.10.7 | Apache-2.0 | 2026-11-12 | `docs/TEST_COVERAGE.md` |
| ruff | 0.14.4 | MIT | 2026-11-13 | One lint pass (below) |
| ViennaRNA | 2.8.0 (wheel `viennarna-2.8.0-cp314-cp314-win_amd64.whl`, SHA-256 `0e68680eff9c1c7865aaa2cfbd313ea3e569e961496f7710e9852651046c252a`) | ViennaRNA Package licence: research, educational and commercial use and modification permitted if the package is not redistributed for a fee and its authors are credited; not redistributed with Synapse | 2026-11-13 | Optional MFE view (`analyses/rna_structure.py`) and the hairpin-screen comparison in `tests/test_rna_secondary_structure.py` (skipped without it). Citation: Lorenz R. et al. (2011) Algorithms Mol Biol 6:26 |

## Static analysis and formatting tools (Master Prompt 7.4)

Evaluated on 2026-11-13:

- **pyright** (in use): type checking at Python 3.10 on every package;
  134 errors at the `ab4bd2a` baseline, 58 now, none new.
- **Ruff** (used once): `ruff check --select F` (the pyflakes rules) over
  the application and the tests found no undefined name and 41 unused or
  duplicated imports (three of them added by the P2-2 move tool); 40 were
  removed with `--fix`, the last disappeared with them; "All checks
  passed". Rendered screens and the suite were unchanged. Not adopted as a
  permanent gate: the project has no Ruff configuration yet, and its
  style rules would flag the hand-formatted tables of the older modules.
- **Black**: not adopted. Reformatting every file would rewrite most lines
  (CRLF files, aligned tables) without any change of behavior, which the
  Master Prompt rules out as an uncontrolled rewrite (1.4, 14.2).
- **mypy**: not adopted; pyright already checks the types.
- **pytest**: not adopted; the `unittest` suite runs on Python 3.10 and
  3.14 without any dependency (`pytest` stays in the `dev` extra).

## Command-line tools

| Tool | Version | Location | Licence | Use in Phase 6 |
|---|---|---|---|---|
| NCBI BLAST+ | 2.17.0+ (build 2025-07-01) | `C:\Program Files\NCBI\blast-2.17.0+\bin` | Public domain (United States Government work) | Detected by the local-BLAST environment checks; no local database is installed (`data/blast_db/` is excluded by design), so no local search was run |
| Strawberry Perl | 5.42.3 | system | Artistic/GPL | Required by `update_blastdb.pl`; not run in Phase 6 |

## Reference data

| Dataset | Source | Version / accession | Obtained | Licence / terms | Checksum | Used by |
|---|---|---|---|---|---|---|
| pUC19 cloning vector, complete sequence (2,686 bp, circular) | NCBI E-utilities, `efetch.fcgi?db=nuccore&id=L09137.2&rettype=fasta&retmode=text` | L09137.2 | 2026-10-30 (audit) and again 2026-11-12; both files identical | NCBI places no restrictions on the use of GenBank data; submitters may hold rights in some records | SHA-256 `656a4ed4597020eb58266ddd4d65f1778abff1a74655abf25ec943dcc14c152a` (LF line ends) | `tests/fixtures/pUC19_L09137.2.fasta`; `PUC19FixtureTests`, `PUC19RestrictionReferenceTests` (132 digests vs Biopython, 0 differences); audit §3.2 |
| REBASE-derived enzyme data bundled with Biopython | Biopython 1.88 `Bio.Restriction` | release as bundled; the API does not expose the REBASE release number | with Biopython | Biopython licence; REBASE data are freely available | – | Enzyme definitions 22/22, Type IIS 4/4, digests, end compatibility |
| REBASE enzyme patterns for EMBOSS | <http://rebase.neb.com/rebase/link_emboss_e> | REBASE version 611 (`emboss_e.611`, dated Oct 29 2026) | 2026-11-13 | Copyright Dr. Richard J. Roberts, all rights reserved; used for comparison only, not stored in the repository | SHA-256 `1253823900d5d7f27383f9f53edaa4330a10a8df4c2bd484791b257cbb57d4a2` (32,393 bytes) | One-off comparison of all 26 enzyme definitions, 0 differences (`docs/METHODS.md`) |
| Primer3 manual worked example | Primer3 manual (CGTGACGTGACGGACT, 50 mM Na+, 50 nM) | as printed in the manual | Previous cycle | – | – | `tests/test_primer_thermodynamics.py` |

No MCP research provider and no project-local skill workflow was executed
in the 2026-11-12 run; the NCBI request above and the python.org downloads
were the only network requests, and they sent no project data. On
2026-11-12 and 2026-11-13 the following further requests were made under
the same policy, again without project data: one live literature search
("CRISPR base editing") through `core.api`, two PubMed metadata look-ups
through the connected PubMed tool (PMID 30587106, 33449100), the REBASE
download above, and one remote NCBI BLAST search of the public pUC19
fragment L09137.2:401-700 against `nt` (plus a first attempt that ended
on a DNS error, see D25).

## Project data and fixtures

| File | Content | SHA-256 |
|---|---|---|
| `data/restriction_enzymes.py` | 22 palindromic enzymes (`synapse-2026-10-29`) and 4 Type IIS enzymes (`synapse-2026-11-11`) | `f95e79b849457aa643816948bb0e9d414f424de1e37700a9a33af28eab6a030f` |
| `examples/dna.fasta` | Example input (`TEST_8_REFERENCE`); origin not recorded | `9060dff9bc85c59405a889bc6262ab69bebb5703acb17510680ca5f430147f9e` |
| `examples/crispr_batch_test.fasta` | CRISPR batch example (`Batch_Test_01` ...) | `c0dd32009fa56f17086d10ae2a73799ee58485ffa3ef79191006071269549f2f` |
| `tests/fixtures/standard_code_reference.py` | Independent copy of NCBI genetic code table 1 | `a11665188fdd262dc1b1abedc1ec89bf34cde8e2544131eec798bf72e2ac8c48` |
| `tests/fixtures/synthetic_circular.gb` | Synthetic 120-bp circular GenBank record ("not a real sequence") | `5e8a7facfbdf7d7d0b0757a7c5c7a97f6733f487412cc325c251e47fcea4101c` |
| `tests/golden/*.json` | 13 golden-master snapshots (internal regression data, not external references) | listed by `sha256sum tests/golden/*.json` |

The checksums of files in Git change when the files change; the Git
history is their authoritative record. They are given here as of
2026-11-12 so that a copy can be checked without Git.
