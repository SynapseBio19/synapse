# Licensing and governance

Prepared in Phase 7 (2026-11-17) for an owner decision. **Decided on
2026-11-20 by the owner: GNU Affero General Public License, version 3 or
later (`AGPL-3.0-or-later`); copyright holder Ali Onur Bayraktar,
privately developed** (section 6). This is an engineering summary, not
legal advice.

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
2. **Copyright holder.** Settled on 2026-11-20: Ali Onur Bayraktar, who
   developed Synapse privately (no university or employer holds rights).
   Copyright line: "Copyright (C) 2026 Ali Onur Bayraktar".
3. **Citation.** Whether to add a `CITATION.cff` so the software can be
   cited (recommended for research software).
4. **Two unused legacy files.** Settled on 2026-11-18 (Phase 8, owner
   mandate).
   - `blast_test.json` is deleted; it remains in the Git history.
   - `TEST_8_REFERENCE_blast.json` is test evidence (a test reads it) and
     moved to `tests/fixtures/legacy_blast/`.
5. **Two unreachable report printers.** Settled on 2026-11-18:
   `print_protein_report` and `print_primer_report` are deleted (no
   caller; they remain in the Git history).

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

## 6. The decision (2026-11-20)

- **Licence:** `AGPL-3.0-or-later`, full text in `LICENSE` (the Free
  Software Foundation's text, SHA-256
  `8486a10c4393cee1c25392769ddd3b2d6c242d6ec7928e1414efff7dfb2f07ef`,
  taken from GitHub's licence API because gnu.org refused automated
  downloads). Declared in `pyproject.toml` (licence and classifier),
  `CITATION.cff` and the README.
- **Owner's goal:** Synapse may become a business later, for example a web
  version on synapse-gen.com. The AGPL lets everyone use, study and change
  Synapse, but whoever distributes a changed version or runs it as a
  network service must publish the source under the same licence. Nobody
  can therefore turn it into a closed product or a closed competing web
  service.
- **The owner keeps every option.** As the only copyright holder, the owner
  can also license Synapse on other terms (for example a paid commercial
  licence) and can choose a different licence for future versions. Copies
  already released under the AGPL stay under the AGPL.
- **Outside contributions.** The owner's freedom to license on other terms
  covers only code the owner holds the rights to. Before merging
  contributions from others, ask the contributor to agree in writing (for
  example a contributor licence agreement, or a statement in the pull
  request) that the owner may also license their contribution on other
  terms. Without that, a contribution can only be used under the AGPL.
- **Compatibility of the material in section 2:** Biopython (BSD-style),
  the bcftools test vectors (MIT/Expat) and BLAST+ (public domain) are
  compatible. primer3-py (GPL-2.0-or-later) is used only by tests and may be
  combined with AGPL-3.0 code under GPL-3.0. ViennaRNA stays optional and is
  not distributed with Synapse.
- **Not done, not required:** licence headers in every source file (the
  `LICENSE` file and `pyproject.toml` declare the licence for the whole
  project; SPDX headers can be added later).
- This summary is not legal advice; for concrete commercial plans a
  lawyer for IT law should review the setup.
