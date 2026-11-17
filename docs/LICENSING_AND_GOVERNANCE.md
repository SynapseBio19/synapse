# Licensing and governance

Prepared in Phase 7 (2026-11-17) for an owner decision. Nothing here
selects a licence: `pyproject.toml` still says "Proprietary - all rights
reserved (decision pending)", and no LICENSE file is added until the owner
chooses one. This is an engineering summary, not legal advice.

## 1. Why the project needs a licence

- Without a licence, copyright law reserves every right to the author.
  Others may read the code if it is published, but they may not copy,
  modify, run in their own work or redistribute it. That rules out
  outside contributors, use in other labs or courses, and publication as
  research software.
- Journals and funders that ask for software availability usually expect
  an OSI-approved licence and a citation file.
- The choice also governs what the project may itself include. The
  dependencies below are compatible with every option in section 3 as
  long as they stay optional.

## 2. Third-party material in or around the repository

| Material | Where | Licence | Consequence |
|---|---|---|---|
| Runtime dependencies | – | none (standard library only) | No constraint on the choice |
| Biopython 1.88 | optional `reference` extra, tests only | Biopython License Agreement (BSD-style) | Permissive; no constraint |
| primer3-py 2.4.0 | optional `reference` extra, tests only | GPL-2.0-or-later | Not imported by the application and not distributed with it. If Synapse ever imported it at run time and was distributed with it, GPL terms would apply to that combination |
| ViennaRNA 2.8.0 | optional `structure` extra, loaded only if installed | ViennaRNA licence: use and modification allowed, no redistribution for a fee, credit required, commercial inclusion by agreement | Keep it optional and not bundled; cite Lorenz et al. 2011 (done in METHODS) |
| bcftools test vectors | `tests/fixtures/bcftools_norm/` (copied) | MIT/Expat (bcftools dual licence) | Keep `SOURCE.md` with the copyright notice; compatible with every option |
| ORFfinder result tables, ExPASy probe measurements, Mutalyzer answers | test fixtures (Synapse's own records of service outputs) | Outputs of public services, recorded with source and date | Kept as test data with provenance; no service code is copied |
| BLAST+ | optional external program, not in the repository | US government work (public domain), NCBI | No constraint |

## 3. Options

| Family | Examples | Others may | Others must | Fits when |
|---|---|---|---|---|
| Permissive | MIT, BSD-3-Clause, Apache-2.0 | use, change and redistribute, also in closed or commercial products | keep the copyright and licence notice (Apache-2.0 also states changes and grants a patent licence) | the goal is the widest use in research and teaching, and closed derivatives are acceptable |
| Weak copyleft | LGPL-3.0, MPL-2.0 | use the code in any project | publish changes to Synapse's own files | changes to Synapse itself should stay open, use inside other software should stay free |
| Strong copyleft | GPL-3.0-or-later | use, change, redistribute | release derived works under the GPL | every derived work, including a future GUI or web front end, should stay open source |
| Network copyleft | AGPL-3.0 | as GPL | also offer the source to users of a hosted service | a future web version (P3-3) must not be run closed by others |
| Source-available, non-commercial | PolyForm Noncommercial | use for research and teaching | not use it commercially | commercial rights are to be kept; not OSI-approved, so many journals and repositories will not count it as open source |
| No licence (current) | "All rights reserved" | nothing beyond reading | – | the code is to stay private |

Practical notes:

- Permissive licences make adoption easiest. MIT is the shortest;
  Apache-2.0 adds an explicit patent grant, which some institutions
  prefer.
- GPL-3.0 keeps derivatives open. It also means that anyone who
  distributes Synapse together with GPL tools, such as primer3-py, faces
  no conflict.
- Changing from a permissive to a copyleft licence later is hard once
  others hold copies. Moving from closed to open is easy.

## 4. Decisions the owner needs to make

1. **Licence.** Which of the families above, and which licence in it.
   After the choice: add `LICENSE`, set `license` in `pyproject.toml` to
   the SPDX identifier, and add a short licence section to the README.
2. **Copyright holder.** `pyproject.toml` names "Onur Bayraktar" as
   author. Confirm the name and year for the copyright line, or name an
   institution if one holds the rights.
3. **Citation.** Whether to add a `CITATION.cff` so the software can be
   cited (recommended for research software).
4. **Two unused legacy files.**
   - `blast_results/TEST_8_REFERENCE_blast.json` (already pending).
   - `blast_test.json` in the project root: a BLAST run from 2026-10-06,
     part of the baseline commit, referenced by no code or document.
   Keep them as historical records, or delete them.
5. **Two unreachable report printers.** `print_protein_report`
   (`ui/protein_menu.py`, about 600 lines) and `print_primer_report`
   (`ui/primer_menu.py`) have no caller. They are older versions of
   reports the menus now draw differently. Keep them, or delete them.

## 5. Governance in place

- **Versioning.** `pyproject.toml`: `0.2.0.dev0`. Recommended: semantic
  versioning, staying below 1.0 until the open items of the Phase 6 and
  Phase 7 reports are settled. Tag a release `v0.2.0` when the licence is
  chosen.
- **Change log.** `docs/CHANGES.md` records every change by phase and
  package, with before, after, reason and evidence. For releases, add a
  short "Release vX.Y.Z" section that points to the packages it contains.
- **Contribution guide.** `CONTRIBUTING.md` lists the setup, the checks
  and the project rules.
- **Repository hygiene.** No generated file is tracked (checked
  2026-11-17). `.gitignore` covers caches, build artifacts, coverage data,
  logs, exports, local BLAST databases, the project-local interpreter and
  the local skills.
- **Continuous integration.** `.github/workflows/ci.yml`. It has not run
  yet because the repository has no remote.
