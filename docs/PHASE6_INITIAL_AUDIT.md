# Phase 6 – independent initial audit

Audit period: 2026-10-30. Scope: the complete project copy
`Project Synapse - Kopie` at commit `0dad3bb` on branch `main`.

This document records what was inspected, what was executed, what was
compared against external reference implementations, and what is defective
or incomplete. It is an audit record, not a change log; no application
source file was modified while it was produced.

**Evidence levels used throughout:**

| Level | Meaning |
|---|---|
| **CONFIRMED** | Reproduced by executing the code in this environment. The observed output is quoted. |
| **SUSPECTED** | Identified by reading the code; a failing input is plausible but was not executed. |
| **UNVERIFIED CONCERN** | A design or scientific question that needs a decision or an external reference, not a demonstrable bug. |
| **NOT A DEFECT** | Examined because it looked wrong, found correct. Recorded so the question is not reopened. |

**Model attribution (required by the Phase 6 brief):** Stage A (baseline,
environment, external reference comparisons for enzymes, digests, Tm,
genetic codes, protein properties, ORF edge cases) was performed with
Opus 5.5. Stages B–D (CRISPR, BLAST, variant, mutation, literature,
remaining modules, cross-module integration, UI and language inventory,
and this document) were performed with Opus 5 after the session was
switched following an interrupted response. No audit task was skipped
because of the switch; §11 lists the one task that was interrupted.
The Git verification, the two repository-housekeeping commits and the
corrections to §1, §10 and §11 were performed with Opus 5.5 after the
model was switched back.

---

## 1. Environment and preservation constraints

| Item | Finding |
|---|---|
| Project root | `C:\Users\bayra\OneDrive\Desktop\Project Snyapse - Kopie` – confirmed the copy, not the original. |
| Branch / HEAD | `main` / `0dad3bb85b597539bb1db067348129738332b0cd` ("Add docs/AUDIT_STATUS.md"). |
| History | 9 commits, read directly from `.git/logs/HEAD`. Stage 0 baseline through Stage 5 all present. Confirmed with `git log` after Git was installed. |
| Working tree | **Not clean at audit start.** `pyproject.toml` carried an uncommitted previous-cycle change and `.claude/` was untracked. An earlier version of this document said "clean"; that was an unverified assumption and is withdrawn. See "Repository state: assumed versus verified" below. |
| **Git executable** | Not available on PATH while the audit was carried out; `.git` was read directly, so no `git diff` review or commit was possible during the audit. **Installed afterwards, on 2026-10-30 (Git for Windows 2.56.0.windows.4), and verified before any commit was made.** |
| Python | 3.14.6 in the project-local `.venv`. `pyproject.toml` declares `requires-python = ">=3.10"`. |
| Python 3.10 | **Not installed here.** The declared minimum version is therefore unverified in this environment. Stage 0 recorded a real 3.10 incompatibility in `analyses/protein.py` (multi-line f-strings), so this is a live risk, not a formality. |
| Installed packages | `pip`, `pyright 1.1.414`, `nodeenv`, `nodejs-wheel-binaries`, `typing_extensions`, plus `biopython 1.88`, `primer3-py 2.4.0`, `numpy 2.5.4` installed for this audit. |
| Dependency policy | Biopython and primer3-py are **reference/test-only**. `pyproject.toml` already declares them under `[project.optional-dependencies] reference`; `dependencies = []` is unchanged. The scientific core remains standard-library only. |
| **BLAST+** | **Installed: `C:\Program Files\NCBI\blast-2.17.0+\bin\` (blastn, blastp, makeblastdb, blastdbcmd).** The previous cycle assumed local BLAST was untestable. It is testable – see §6. |
| `data/blast_db/` | Absent by design. Six catalog entries (`nt`, `nr`, `refseq_*`, `swissprot`) report `NOT_INSTALLED` cleanly. No large download was performed or is required. |
| `data/local_test/` | **A real 3-sequence BLAST nucleotide database is bundled** (`synapse_test.n*`) plus `query.fasta`. **No test uses it.** |
| bcftools / tabix / samtools / primer3_core | None on PATH. `bcftools` has no supported native Windows build; a `bcftools norm` comparison needs WSL or a Linux host. |

### Repository state: assumed versus verified

Two sources describe the repository at the start of the audit. They
disagree, and only the second is evidence.

| | Assumed before Git was available | Verified with Git 2.56.0.windows.4 on 2026-10-30 |
|---|---|---|
| Source | Status reported by the session environment at start-up. No Git binary existed, so it could not have been computed. | `git log`, `git status`, `git diff`, `git stash list`, run once Git for Windows was installed. |
| HEAD | `0dad3bb`, read from `.git/refs/heads/main` | `0dad3bb` — confirmed |
| Tracked files | Reported "clean" | **`pyproject.toml` modified** (+10/−1): adds `pyright[nodejs]==1.1.414` to the `dev` extras and a `[tool.pyright]` table targeting Python 3.10. Last written 2026-10-30 02:11:17, before this session began, and already present at this session's first read of the file. Previous-cycle work, not Phase 6 work. Runtime dependencies (`dependencies = []`) and the `reference` extras unchanged. |
| Untracked files | Not reported | **`.claude/`** — 101 files (1.21 MB), consisting only of the eleven K-Dense skills under `.claude/skills/`. Present at this session's first file listing. Plus the two Phase 6 documents. |
| Staged changes / stash / deleted files | Not reported | None / none / none |
| Branches and tags | Not reported | Branch `main` only; tag `baseline-audit` |

Neither difference was introduced by this audit, and neither affects any
finding in this document: `pyproject.toml` contains no application logic,
and `.claude/skills/` lies outside every path that the application, the
test discovery (`tests/`) and pyright's `include` list read.

Resolution, 2026-10-30 (approved repository housekeeping; no application
source file touched):

- `5242fcb` commits the pre-existing `pyproject.toml` change, unchanged and on its own.
- `5d8c73a` adds `/.claude/skills/` to `.gitignore`. All eleven skills remain on disk; the anchored rule leaves other `.claude/` content visible to Git.
- This document and `docs/FEATURE_STATUS.md` are then committed together.

### Project structure actually present

```
main.py                 633 lines   terminal entry point, module dispatch
analyses/   12 modules   36,096 lines total (blast 4,650; primer 5,412;
                                     cloning 4,351; mutation 4,184;
                                     protein 3,777; rna 2,617; crispr 2,574;
                                     orf 2,587; variant_analysis 2,322;
                                     restriction 2,230; sequence 1,768;
                                     literature 1,624)
core/        8 modules    1,806 lines   api, pipeline, models, report,
                                        seqtools, orf_context, validation,
                                        logging_setup
fileio/      3 modules      255 lines   fasta, genbank, loader
ui/          2 modules      772 lines   menu, report
data/        3 modules      400 lines   restriction_enzymes, amino_acids,
                                        genetic_code
tests/      11 files        ~1,200 lines  109 tests + 12 golden scenarios
docs/        3 files        AUDIT_STATUS, CHANGES, METHODS
LEGACY/                  pre-project scripts, excluded from packaging
Synapse_Website/         static placeholder page, out of Phase 6 scope
blast_results/           one stored BLAST run (see D7)
examples/                dna.fasta (14 records), crispr_batch_test.fasta
```

Approximately **1,100 functions and classes** were inventoried across the
twelve analysis modules. Every module was opened; every module's public
entry points and menu dispatcher were read; all twelve module self-test
suites and the full unit-test suite were executed.

---

## 2. Baseline test run

```
python -m unittest discover -s tests -t .
Ran 109 tests in 16.591s
OK
```

All twelve built-in module self-tests also pass when invoked directly:

| Module self-test | Result |
|---|---|
| `crispr.run_self_test()` | 28/28 pass |
| `blast.run_self_test()` | 33/33 pass |
| `variant_analysis.run_variant_analysis_self_test()` | 27/27 pass |
| `mutation.run_mutation_self_test()` | 44/44 pass |
| `sequence.run_sequence_self_test()` | all pass |
| `literature.self_test()` | all pass (offline fixtures) |
| `restriction`, `cloning`, `primer`, `protein`, `rna`, `orf` | pass via `tests/test_builtin_selftests.py` |

**A passing suite is not evidence of scientific correctness.** Every defect
in §4 and §5 below exists in code that the suite reports as green. The
suite's blind spots are quantified in §8.

---

## 3. External reference comparisons performed in this audit

These are genuine comparisons against independent implementations, executed
in this environment. They are distinct from internal unit tests and from
golden-master snapshots.

### 3.1 Restriction enzyme definitions – CLOSES the open P0-4 item

All 22 entries in `data/restriction_enzymes.py` were compared against
`Bio.Restriction` (REBASE-derived EMBOSS data bundled with Biopython 1.88)
on recognition sequence, top-strand cut offset, bottom-strand cut offset,
overhang polarity and overhang length.

**Result: 22/22 agree on all five properties.** Including the two entries
the previous cycle corrected or flagged: `ClaI` `AT^CGAT` with a 2 nt 5′ CG
overhang, and `NdeI` `CA^TA_TG` with a 2 nt 5′ TA overhang.

### 3.2 In-silico digests on a real plasmid – CLOSES the open P1-4 item

Reference: **pUC19, NCBI accession L09137.2**, 2,686 bp, circular,
downloaded 2026-10-30 from NCBI E-utilities;
FASTA SHA-256 `656A4ED4597020EB58266DDD4D65F1778ABFF1A74655ABF25EC943DCC14C152A`.
Kept outside the repository (scratchpad only).

`analyses.restriction.perform_digest` was compared with
`Bio.Restriction.<enzyme>.search()` for all 22 enzymes in both linear and
circular topology, on three sequence rotations:

| Sequence | Comparisons | Differences |
|---|---|---|
| pUC19 as published | 44 | **0** |
| pUC19 rotated so an EcoRI site spans the origin | 44 | **0** |
| pUC19 rotated so a BamHI cut falls exactly at the origin | 44 | **0** |

Cut coordinates and the full sorted fragment-length lists match in every
case. The `wraps_origin` handling introduced in Stage 3 is correct.

### 3.3 End compatibility

The full 22 × 22 = 484-pair compatibility matrix from
`analyses.cloning.ends_compatible` was compared against Biopython's
`ovhg`/`ovhgseq` criterion (equal overhang polarity and equal overhang
sequence). **0 mismatches.** P0-3 confirmed.

### 3.4 Melting temperature – CLOSES the open P1-3 item

`analyses.primer.calculate_tm_nearest_neighbor` vs
`primer3.calc_tm(tm_method="santalucia", salt_corrections_method="santalucia")`:
8 oligos × 5 buffer conditions (Na⁺ 50–100 mM, Mg²⁺ 0–5 mM, dNTP 0–0.8 mM,
primer 50–250 nM) = **40 comparisons, maximum |Δ| = 0.003 °C.**

This satisfies the original acceptance criterion ("±1 °C against
Primer3/Biopython at three Mg levels") with a margin of more than two
orders of magnitude. The residual 0.003 °C is the gas-constant rounding
already documented in `docs/METHODS.md`.

### 3.5 Genetic codes

`core.seqtools.GENETIC_CODES` vs `Bio.Data.CodonTable` for NCBI tables
1, 2, 4 and 11: all 64 codon assignments identical in every table, and the
start-codon sets identical (table 1: ATG/CTG/TTG; table 2: ATA/ATC/ATG/ATT/GTG;
table 4: 8 starts; table 11: 7 starts).

### 3.6 Protein physicochemical properties

`analyses.protein` vs `Bio.SeqUtils.ProtParam.ProteinAnalysis`, on ubiquitin
(75 aa) and hen lysozyme (129 aa):

| Property | Ubiquitin: Synapse / Biopython | Lysozyme: Synapse / Biopython |
|---|---|---|
| Isoelectric point | 6.5616 / 6.56 | 9.3238 / 9.32 |
| GRAVY | −0.48947 / −0.489 | −0.47209 / −0.472 |
| Instability index | 36.055 / 36.06 | 16.091 / 16.09 |
| Aromaticity | 0.039474 / 0.0395 | 0.093023 / 0.093 |
| Extinction (reduced / oxidised) | 1490 / 1490 | 37470, 37970 / 37470, 37970 |
| Molecular weight | 8564.84 | 8564.74 (Biopython) |

The molecular-weight difference is a residue-mass-table difference;
Synapse's 8564.84 matches the ExPASy ProtParam value already cited in
`docs/METHODS.md`. All other values agree to the precision Biopython prints.

### 3.7 Variant left-alignment

`analyses.variant_analysis.left_normalize_variant` on constructed repeat
contexts, checked against the Tan et al. (2015) algorithm by hand:

| Reference | Input | Output | Correct? |
|---|---|---|---|
| `GGGAAAAAGGG` | `7 AA>A` | `3 GA>G` | Yes – shifted to the left edge of the A-run and re-anchored |
| `GGGAAAAAGGG` | `5 A>AA` | `3 G>GA` | Yes |
| `GTCACACACACAG` | `7 CA>C` | `7 CA>C` | Yes – after trimming this is a 1 nt A deletion preceded by C, which cannot shift |
| `GTCACACACACAG` | `2 T>TCA` | `2 T>TCA` | Yes – insertion's last base A ≠ preceding base T |
| `GGGCACACACAGGG` | `8 ACAC>A` | unchanged + `REF_MISMATCH` | Yes – deliberately wrong REF, correctly flagged |

`trim_common_bases` behaves correctly (right-trim then left-trim with
re-anchoring): `(5,"CTAG","CTG") → (7,"AG","G")`, `(1,"AAAT","AAT") → (3,"AT","T")`.

**This is a hand-verified algorithm check, not a `bcftools norm`
comparison.** bcftools remains unavailable (§1). The P0-7 status should be
upgraded from "not independently validated" to "algorithm verified against
the published method on constructed cases; tool comparison still open".

---

## 4. Confirmed defects

### D1 – ORF translation ignores the policy-selected initiator set
**`analyses/orf.py:2651`** · CONFIRMED · scientific

`translate_orf_sequence` hardcodes `ALL_START_CODONS` (ATG/GTG/TTG) as the
initiator set whenever `genetic_code` is 1, instead of the NCBI table-1
start set (ATG/CTG/TTG) or the policy actually used for the search.

With `start_codon_policy="table_starts"` an ORF beginning `CTG` is found,
but:

```
O.translate_orf_sequence(orf["sequence"])      -> "LAAAA"
core.api.analyze_protein_from_cds(orf["sequence"]).data["protein"] -> "MAAAA"
```

Two code paths give different first residues for the same ORF. The ORF
menu's own translation preview disagrees with the protein module and with
the shared API.

### D2 – ORF translation silently drops the ORF's own genetic code
**`analyses/orf.py:2637`** · CONFIRMED · scientific

`translate_orf_sequence(orf_sequence, initiator_start=True, genetic_code=1)`
defaults to table 1. ORF dictionaries carry `genetic_code`, but the
function does not read it, so any caller that omits the argument
mistranslates a non-standard-table ORF. For a table-2 ORF
(`ATG` + `TGA`×25 + `TAG`):

```
translate_orf_sequence(seq, genetic_code=2) -> "MWWWWWWW..."   (TGA = Trp)
translate_orf_sequence(seq)                 -> "M"             (TGA = stop)
```

The protein is truncated to a single residue with no warning.
`analyze_selected_orf` does pass the table; the hazard is the default.

### D3 – Restriction map and site finder have no topology parameter
**`analyses/restriction.py:201, 469`** · CONFIRMED · scientific

`find_restriction_sites(sequence, enzyme_name)` and
`create_restriction_map(sequence, enzyme_names)` take no `topology`
argument and call `find_restriction_site_matches` with its default
`topology="linear"`. On pUC19 rotated so an EcoRI site spans the origin:

```
perform_digest(rot, ["EcoRI"], topology="circular")  -> site found, correct fragments
create_restriction_map(rot, ["EcoRI"])               -> []
```

`analyze_restriction(..., topology="circular")` builds its
`restriction_map` through `create_restriction_map`, so the circular digest
is right while the map, the site table shown to the user, and every
consumer of the map are wrong for origin-spanning sites.

### D4 – Cloning module repeats the P0-5 ambiguity defect
**`analyses/cloning.py:276`** · CONFIRMED · scientific

`cloning.find_restriction_sites` matches with IUPAC tolerance and returns a
bare position list with **no CONFIRMED/POSSIBLE distinction**:

```
cloning.find_restriction_sites("A"*10 + "N"*6 + "A"*10, "GAATTC") -> [11]
```

A run of N produces a phantom EcoRI site. The identical defect was fixed in
`restriction.py` as P0-5 (`include_possible=False` by default) but
`cloning.py` keeps its own independent implementation, so
`analyze_enzyme`, `build_restriction_map`, `analyze_enzyme_pair`,
`automatic_cloning_planner` and every strategy score built on them can be
based on sites that do not exist.

### D5 – Cloning module has no circular-topology support
**`analyses/cloning.py`** · CONFIRMED · scientific

Neither `find_restriction_sites`, `analyze_enzyme`, `build_restriction_map`
nor `_digest_fragments`' callers accept or propagate a topology for site
discovery (`_digest_fragments` itself has a `topology` parameter, but the
sites fed to it are found linearly). A circular vector — the normal case in
cloning — is analysed as if it were linear, so sites and junctions spanning
the origin are invisible to the planner.

### D6 – Programmatic local BLAST is unreachable
**`analyses/blast.py:3379`** · CONFIRMED · functional

```python
db_path = Path(local_database_path).expanduser().resolve()
if not db_path.exists():
    raise BlastAnalysisError(f"Local BLAST database path not found: {db_path}")
```

A BLAST database is addressed by a **path prefix**, not a file. For the
bundled database:

```
Path(".../data/local_test/synapse_test").exists()      -> False
Path(".../data/local_test/synapse_test.nin").exists()  -> True
run_blast_analysis(..., engine="local",
                   local_database_path=".../synapse_test")
   -> BlastAnalysisError: Local BLAST database path not found
```

So `run_blast_analysis(engine="local", local_database_path=…)` — the
documented programmatic entry point, and the one a future UI or script
would call — can never succeed with a valid database. Only the menu path,
which goes through `LocalBlastDatabaseManager.scan_database`, works.
Bypassing the gate proves the rest of the pipeline is sound (§6).

### D7 – A stale, scientifically wrong BLAST artifact is shipped
**`blast_results/TEST_8_REFERENCE_blast.json`** · CONFIRMED · scientific

The stored run still carries pre-P0-1 values:

```
query.length                      : null
hits[0].best_coverage_percent     : 100.0        (32 bp HSP on a 99 bp query)
hits[0].hsps[0].evalue            : 7.45249
hits[0].hsps[0].query_from/to     : null
hits[0].evidence                  : absent
```

It also uses an older HSP key schema (`query_from`/`query_to`) than the
current `BlastHSP` (`query_start`/`query_end`/`subject_start`/`subject_end`),
i.e. the file cannot be re-read by today's code. **No test references it.**
`docs/AUDIT_STATUS.md` cites this file as the evidence for the P0-1 fix
while the file itself still shows the defective output.

### D8 – CRISPR reports no cut site
**`analyses/crispr.py`** · CONFIRMED · scientific / functional

`GuideCandidate` has 28 fields and none of them is a cut or cleavage
position. Across the whole module `cut_site` and `cut_position` occur **0**
times; `guide_to_dict` and the CSV/JSON exports contain no cut key. The
only cleavage information is the free-text `CRISPRSystem.cut_model` string.
P2-3 item 1 ("correct CRISPR cut-position handling") is not implemented,
and the menu offers "Target Discovery" and "Guide Analysis & Ranking"
without ever stating where the enzyme is modelled to cut.

### D9 – The on-target site is scored as an off-target
**`analyses/crispr.py:1146, 846`** · CONFIRMED · scientific

The local specificity screen does not exclude the guide's own intended
site. Designing against a sequence and then supplying that same sequence as
the local reference — exactly what the menu's "Off-Target Analysis" invites —
gives:

```
specificity_profile: exact_matches = 1, total_hits = 1
_component_specificity -> 60.0      (100.0 − 1 × 40.0)
```

A perfectly specific guide is penalised 40 points for matching its target,
and that penalty enters `rank_guides` with weight 0.10. Guide ranking is
therefore distorted in the one workflow the module advertises for
specificity. P2-3 item 3 not implemented.

### D10 – Five CRISPR configuration fields are dead
**`analyses/crispr.py:261`** · CONFIRMED · functional

`include_low_quality`, `rank_by_gc`, `editing_window_start`,
`editing_window_end` and `base_editor` appear in `DesignConfig` but are
never read anywhere in the module (0 attribute references). Behavioural
proof:

```
include_low_quality=True  -> candidates=3, ids G0001..G0003
include_low_quality=False -> candidates=3, ids G0001..G0003   (identical)
rank_by_gc=True / False   -> scores [79.239, 75.788, 8.042]   (identical)
```

`base_editing_window()` returns `window or definition["window"]`, so a
user-configured editing window in `DesignConfig` is ignored in favour of
the hardcoded `BASE_EDITOR_DEFINITIONS` window. These are advertised
options that do nothing.

### D11 – CRISPR exports default outside the project
**`analyses/crispr.py:1668`** · CONFIRMED · architectural

```
_default_export_directory() -> C:\Users\bayra\OneDrive\Desktop
```

Guide exports land on the user's Desktop (here, inside a cloud-synced
folder) rather than in a project-controlled export directory. `.gitignore`
already reserves `/exports/`, which is what the report exporter uses.

### D12 – Variant normalisation fails silently on contig-name mismatch
**`analyses/variant_analysis.py:979, 958`** · CONFIRMED · scientific

`check_reference_alleles` looks the contig up by exact name
(`contig_references.get(variant.chromosome)`), so the ubiquitous `chr1`
vs `1` difference is unhandled in both directions:

```
reference key 'chr1', variant 'chr1' -> checked=1, mismatched=0
reference key 'chr1', variant '1'    -> checked=0, mismatched=0, not_checked=1
reference key '1',    variant 'chr1' -> checked=0, mismatched=0, not_checked=1
```

Worse, `normalize_variants` then returns the variant **unchanged with no
flag**:

```
normalize_variants([Variant("1", 7, "AA", "A")], {"chr1": "GGGAAAAAGGG"})
   -> 7 AA>A, validation_flags: None        (expected 3 GA>G)
```

while `core/api.analyze_vcf` records `parameters["normalized"] = True`
because a reference mapping was supplied. A user with a `1`-style reference
and a `chr1`-style VCF is told the variants were normalised and REF-checked
when neither happened. `mismatched: 0` reads as success; the real signal is
buried in `not_checked`, which no caller surfaces.

### D13 – Two live `normalize_dna` functions with different contracts
**`core/validation.py:1` vs `core/seqtools.py:36`** · CONFIRMED · architectural

```
core.validation.normalize_dna("aU g") -> "AUG"   (keeps U)
core.seqtools.normalize_dna("aU g")   -> "ATG"   (maps U -> T)
core.validation.validate_dna("ACGU")  -> False
```

Same name, different semantics, both in `core/`, both in use: `main.py`
imports the `validation` pair, `core/api.py` and the newer modules import
the `seqtools` one. An RNA FASTA is rejected by the menu but silently
accepted as DNA through the shared API. `core/seqtools.py`'s own docstring
says the older modules "will be migrated without behavior change"; this
divergence means the two contracts are not interchangeable.

### D14 – No alphabet detection: a peptide loads as DNA
**`fileio/`, `core/validation.py`** · CONFIRMED · scientific

Every letter of `MKVHDANWYSK` is a valid IUPAC *nucleotide* code, so:

```
load_sequences(prot.fasta)  -> 1 record, sequence "MKVHDANWYSK", no note
core.validation.validate_dna("MKVHDANWYSK") -> True
```

A protein FASTA is accepted as a DNA record and flows into GC content, ORF
search, restriction analysis and primer design without a single warning.
`SequenceRecord` has an `alphabet` field defaulting to `"DNA"` that nothing
ever sets or checks.

### D15 – An advertised main-menu feature is a placeholder
**`main.py:762`** · CONFIRMED · functional

`[11] AI Experiment Planner` prints "AI Experiment Planner is coming soon."
It is the only main-menu entry with no implementation.
`docs/AUDIT_STATUS.md` lists it as P3-1 DEFERRED, but the menu advertises
it as a feature.

### D16 – Duplicate and shadowed primer functions
**`analyses/primer.py:1782, 3699`** · CONFIRMED · architectural

`calculate_tm_nearest_neighbor` is defined twice at module level. The
second definition is a backward-compatibility adapter that captures the
first as `_calculate_tm_nearest_neighbor_v3_raw` and replaces the public
name — deliberate, documented, and correct, but it means the name a reader
finds first is not the name that is called. The live `_v5`–`_v11` ranking
wrapper chain (eight generations of `*_pair_quality` and
`*_rank_professional_primer_pairs`) also remains, as `docs/AUDIT_STATUS.md`
P2-1 records. `analyses/primer.py` is 5,412 lines, the largest module in
the project, and its public surface no longer reflects its structure.

### D17 – The terminal interfaces share no visual language
**all modules** · CONFIRMED · interface

Mechanically counted separator literals (`"=" * n`, `"-" * n`):

| Width | Occurrences | Width | Occurrences |
|---|---|---|---|
| 42 | 5 | 70 | 130 |
| 45 | 108 | 72 | 35 |
| 58 | 1 | 78 | 34 |
| 60 | 51 | 88 | 5 |
| 65 | 2 | 90 | 5 |
| 68 | 11 | | |

**Eleven distinct widths.** Per module: BLAST 90, literature 88, CRISPR and
restriction 78, primer 72, mutation / variant / rna / protein / orf 70 with
45 for sub-sections, cloning 60 + 70, `main.py` 60 **and** 68, `ui/menu.py` 60.
Two modules (`restriction`, `orf`) mix three or more widths internally.

Navigation also diverges: `[0]` means **Back** in all eleven sub-module
menus but **Exit** in the main menu (`ui/menu.py`); `blast` and `crispr`
use `[0] Cancel` in sub-dialogs and `[0] Back` in menus; letter keys
(`[P]`, `[R]`, `[S]`, `[F]`) exist only in the main menu and the sequence
selector.

### D18 – Exported results still carry the "confidence" vocabulary
**`analyses/orf.py:810`, exports** · CONFIRMED · editorial / scientific

P2-5 renamed the displayed ORF label from "Confidence" to "Heuristic
score", but the data model and every export still use the old vocabulary:

```
orf["confidence"] = {"score": 75, "classification": "Moderate"}
```

`confidence` occurs 49 times in `analyses/`. JSON and Markdown reports,
`compare_orf_candidates` (`confidence_score`, `confidence_classification`)
and the golden snapshots all present a heuristic ranking aid under a name
that implies a probability. The display was fixed; the interface was not.

### D19 – The pipeline's ORF runner cannot reproduce the menu's search
**`core/pipeline.py:627`** · CONFIRMED · architectural

`SynapsePipeline.run_orf(include_partial=True, min_length=None,
enhanced=True, top_n=5)` exposes none of `start_codon_policy`, `nested`,
`ranking` or `genetic_code`. The ORF menu's "[6] Search Settings" can
change all four, and `core/api.find_orfs` accepts all four, but the
pipeline's own runner is stuck on the defaults — so a programmatic
pipeline run cannot reproduce an interactive session's ORF result.

`PipelineContext.completed` is initialised with six keys while
`store_result` writes `self.completed[result.module] = True` for any
module name, so the dict silently grows to whatever has been stored and
`status()` returns a history-dependent shape, while `print_status` prints a
fixed list of twelve labels. Not a crash; an inconsistent contract.

---

## 5. Suspected defects and unverified concerns

| ID | Location | Level | Finding |
|---|---|---|---|
| S1 | `analyses/orf.py:483` | SUSPECTED | `analyze_kozak_context` is called with the *strand-local* sequence and index, so for a minus-strand ORF the reported Kozak context belongs to the reverse-complement coordinate frame. Observed output is `classification: "Unavailable"` for the minus ORFs tested (the window ran off the sequence end), so no wrong value was produced here, but the coordinate frame is undocumented and a longer flanking region would produce a context whose position cannot be mapped back to the parent sequence. |
| S2 | `analyses/orf.py:942` | SUSPECTED | `_keep_longest_per_stop` keys on `(strand, frame, three_prime, status)`. Including `status` means a complete ORF and a partial ORF sharing a 3′ end are both kept, which is defensible, but the key uses `orf["stop"]` for `+` and `orf["start"]` for `−` as the 3′ end — correct — while `frame` for the minus strand is derived from the reverse-complement offset, so two minus ORFs in the same biological frame but found from different rotations would not collide. Needs a targeted test. |
| S3 | `analyses/crispr.py:1326` | UNVERIFIED CONCERN | Prime-editing PBS/RTT construction is anchored on `mutation.position` and emits a candidate only when the requested lengths fit. The construction is labelled heuristic and carries explicit notes, which is honest, but there is no check that the RTT actually contains the intended edit — `edit_in_rtt` only tests `len(rtt) == rtt_length`. The field name overstates what is verified. |
| S4 | `analyses/crispr.py:135` | UNVERIFIED CONCERN | `DEFAULT_MIN_GC = 30`, `DEFAULT_MAX_GC = 70`, `DEFAULT_MAX_HOMOPOLYMER = 4`, `DEFAULT_SEED_LENGTH = 8`, base-editing window `(4, 8)`. These are plausible community conventions but carry no citation in the module or in `docs/METHODS.md`. They should either be referenced or explicitly labelled as Synapse conventions. |
| S5 | `analyses/crispr.py:642` | UNVERIFIED CONCERN | The poly-T flag fires on any run of ≥4 T anywhere in the guide. The U6-promoter consideration that motivated P2-3 item 2 concerns a TTTT terminator motif; the current flag is a generic homopolymer warning, is position-independent, and the string "U6" appears nowhere in the module. The item is addressed only incidentally. |
| S6 | `analyses/blast.py` | SUSPECTED | `NCBIBlastClient` exposes `submit`/`status`/`fetch_xml`/`run`, i.e. the primitives needed to resume a search from a saved RID, and "rid" occurs 40 times, but "resume" occurs 0 times and no menu entry offers it. P1-7's "RID resume deferred" is accurate; the capability is one small function away from existing. |
| S7 | `analyses/blast.py:3385` | SUSPECTED | On the user-supplied-path branch, `database_details` is populated with `title`, `database_date`, `blastdb_version`, `sequence_count` all `None` because `blastdbcmd -info` is not invoked. Confirmed in the §6 run. Provenance is weaker on that path than on the managed path. |
| S8 | `core/report.py:72` | SUSPECTED | `to_jsonable` catches bare `Exception` around `value.to_dict()` and falls through to `repr(value)`. A payload whose `to_dict` raises is silently serialised as a Python repr string, which is unparseable downstream and gives no warning in the report. |
| S9 | `analyses/mutation.py` | UNVERIFIED CONCERN | The global Needleman–Wunsch alignment is O(n·m) in time *and memory*; the only mitigation is a note above 2,000 nt (P1-10). P2-6 records the 10 kb target as unmet. For a 10 kb × 10 kb comparison the DP matrix is ~10⁸ cells, which will exhaust memory on a typical machine rather than merely being slow. The note should state the limit, not just flag it. |
| S10 | `fileio/genbank.py` | SUSPECTED | Location parsing covers simple spans, `complement()` and `join()`; other forms keep only the raw string. A record using `order()`, `<`/`>` partial bounds or a remote accession reference therefore yields a feature with no usable interval, silently. The limitation is documented in METHODS but not surfaced as a per-record note. |
| S11 | `analyses/literature.py` | UNVERIFIED CONCERN | Four providers (PubMed, Europe PMC, OpenAlex, Crossref) each fail independently, and `_enrich_missing_metadata` fills gaps across sources. Dedup (`deduplicate_articles`) and relevance ranking were exercised only through the module's offline fixtures. No live call was made in this audit, by design. Real-world behaviour on rate limits, partial outages and cross-source duplicate DOIs is therefore unverified. |
| S12 | project-wide | UNVERIFIED CONCERN | `pyproject.toml` sets `[tool.pyright] pythonVersion = "3.10"`, `typeCheckingMode = "standard"`. `pyright` is installed but **was not run** in this audit (read-only scope; it writes caches). Running it is the cheapest available check for the 3.10 compatibility risk in §1. |
| S13 | `main.py` | SUSPECTED | `main.py` is a top-level script with module-level control flow (`sequences = choose_sequence_file()` at import, `while True:` loops at module scope) and no `if __name__ == "__main__":` guard. `tests/test_main_smoke.py` therefore has to drive it by patching `builtins.input`. Any import of `main` starts the application. |
| S14 | `analyses/*` | UNVERIFIED CONCERN | Eleven of twelve analysis modules own their terminal menus (`run_*_ui`, `*_menu`), so `print()`/`input()` and scientific computation live in the same file: 1,966 `print(` and 126 `input(` calls across `analyses/`. `core/api.py` provides a clean non-interactive layer for ten analyses, which is the P2-2 compromise, but the modules themselves are not UI-free and `core/api.py` does not yet cover cloning, variant filtering, mutation reporting, BLAST or literature. |

Found while verifying the D12 correction (2026-10-31); neither is in D12's
scope and neither was changed:

| ID | Location | Level | Finding |
|---|---|---|---|
| S15 | `analyses/variant_analysis.py`, `annotate_from_mutation_analysis` (call at lines 2026–2027 after D12) | SUSPECTED | pyright reports `No parameter named "reference_cds"` and `"variant_cds"` on a call into the mutation module. If that call path is reached it would raise `TypeError`. Present identically at `824ca5f`; not reproduced by execution. |
| S16 | `analyses/variant_analysis.py`, `left_normalize_variant` | CONFIRMED | The normalized copy is built with `annotation=dict(variant.annotation)`, a shallow copy, and `REF_MISMATCH` is then added with `setdefault(...).append`, so the flag is also written into the *input* variant's `validation_flags` list when that list already exists. Reproduced: input flags `["PRE_EXISTING"]` become `["PRE_EXISTING", "REF_MISMATCH"]` after `left_normalize_variant(Variant("chr1", 1, "TTT", "T"), "GGGAAAAAGGG")`; an input without a pre-existing list stays clean. D12's new flags use copy-on-write and do not have this problem; the pre-existing `REF_MISMATCH` path was left unchanged. (Listed here with the S-items because it was found after the audit; it is a confirmed defect.) |

Found during the D4 work (2026-10-31); outside D4's scope and not changed:

| ID | Location | Level | Finding |
|---|---|---|---|
| S17 | `analyses/cloning.py`, `_selected_digest` | CONFIRMED | The cloning digest uses each site's recognition **start** as the cut position instead of the top-strand cut coordinate, so fragment lengths are off by the enzyme's cut offset. Reproduced with EcoRI (`G^AATTC`) + PstI (`CTGCA^G`) on a 162 bp sequence: cloning reports cuts `[51, 107]`, fragments 51/56/55 bp; `restriction.perform_digest` reports cuts `[51, 111]`, fragments 51/60/51 bp. EcoRI agrees only because its cut offset happens to equal 1. |
| S18 | `analyses/cloning.py`, `automatic_cloning_planner` | CONFIRMED | Runtime: on a 1.3 kb vector and a 0.3 kb insert the planner takes 100–110 s with 10 enzymes and well over 2 minutes with all 22 (environment-dependent; measured 135 s on a lightly loaded run), because every strategy re-scans both molecules (about 12 site searches per evaluated strategy, ~63,500 strategies for 22 enzymes). Results are correct; the cost is repeated work. Related to P2-6. |

Found during the D1 work (2026-10-31); outside D1's scope and not changed:

| ID | Location | Level | Finding |
|---|---|---|---|
| S19 | `analyses/orf.py`, `classify_start_codon`, `calculate_cds_quality`, `calculate_orf_confidence` | CONFIRMED | Start-codon *classification* and the heuristic scores still use the fixed ATG/GTG/TTG set rather than the genetic code. A CTG ORF found with `table_starts` gets `classify_start_codon("CTG") == "None"`, no "Alternative start codon" reason, and a CDS-quality score of 60 versus 70 for an otherwise identical GTG ORF; table-2 starts such as ATA/ATT are treated the same way. Translation is not affected (D1). |

### Examined and found correct (not defects)

- `crispr.list_crispr_systems()` returns `CRISPRSystem` objects, not names; all three call sites use them correctly.
- Cas12a 5′-PAM geometry: PAM at 41–44, guide at 45–67, `target_sequence` ordered PAM+guide — verified against the sequence.
- CRISPR coordinate round-trip on both strands: 0 mismatches over all candidates in three synthetic sequences, including a sub-region with `region_start=61` (parent-consistent) and a minus-strand ORF region (5–100, forward region identical to the parent slice). The P0-2 fix holds.
- ORF minus-strand coordinates: `start`/`stop` are 0-based inclusive forward-strand with `start ≤ stop`, and `reverse_complement(dna[start:stop+1])` reproduces `sequence` exactly.
- ORF-change invalidation: storing `rna`, `protein`, `crispr`, `blast` results and then changing the selected ORF clears all four (`{} ` afterwards). P1-2 holds.
- The RNA selected-ORF-versus-whole-sequence rule is intact: with no ORF selected the source is the whole sequence and the label is "Selected DNA"; with an ORF selected the source is the ORF sequence and the label is "Selected ORF"; `rna.resolve_rna_source` additionally reports `mode: "orf"`.
- `fileio` rejects an empty file and a header-less file with a clear message, and content-sniffs rather than trusting the extension.
- Report export rejects any suffix other than `.json`/`.md`.
- Literature date bounds: `2020-03` → lower `2020-03-01`, upper `2020-03-31`; `2020` → upper `2020-12-31`.
- Evidence bands are gated on E ≤ 1e-3 (`DEFAULT_EVIDENCE_EVALUE_THRESHOLD = 0.001`) with a `NOT_SIGNIFICANT` band, and the threshold is parameterised.
- BLAST hit identity of 100 % alongside 98.04 % coverage (§6) is correct BLAST semantics — identity is per-HSP over the aligned length — and matches how NCBI reports it.

---

## 6. Local BLAST: first real execution

Previously recorded as untestable. BLAST+ 2.17.0+ is installed and the
repository bundles a 3-sequence nucleotide database, so the local path was
executed end to end through Synapse's own engine (bypassing D6's faulty
path gate by constructing the `LocalBlastDatabase` record directly):

```
query  : Synapse_Test_Query, 51 nt, DNA
program: blastn      database: data/local_test/synapse_test
engine : Local BLAST+        version: BLASTN 2.17.0+

hits: 2
  Test_Gene_001  coverage 100.00 %  identity 51/51  E 1.19e-25  bits 95.30
  Test_Gene_002  coverage  98.04 %  identity 50/51  E 4.28e-25  bits 93.45

reported_parameters: expect=10, sc-match=1, sc-mismatch=-2,
                     gap-open=0, gap-extend=0, filter=L;m;
evidence: Test_Gene_001 VERY_STRONG_SEQUENCE_EVIDENCE (score 100.0)
```

This validates, on real BLAST output rather than synthetic XML:

- the **P0-1 coverage definition** — 51/51 = 100 %, 50/51 = 98.04 %, computed
  from query intervals against the query length;
- XML parsing, hit and HSP construction, and `BlastRunMetadata` provenance
  (program, database, version, effective parameters, local path, timestamps);
- evidence banding with the significance gate;
- the text report renderer.

Two gaps follow: `database_details` is unpopulated on this path (S7), and
**no automated test exercises any of it** — neither BLAST+ nor the bundled
database is referenced anywhere in `tests/`.

---

## 7. Cross-module integration

| Workflow | Status | Evidence |
|---|---|---|
| DNA → validate → ORF → select → RNA → protein | **Works** | Source selection, labels and invalidation all verified (§5 "found correct"). `tests/test_main_smoke.py` drives it end to end. |
| Selected-ORF genetic code → protein | **Works via the menu, fragile via the API** | `main.py` warns when a non-table-1 ORF is selected; `analyze_complete_protein` receives the table. But D2 makes any direct `translate_orf_sequence` call a silent mistranslation. |
| Selected ORF → CRISPR | **Works** | Region 5–100 on the parent, correct strand, forward region matches the parent slice, 0 coordinate mismatches. |
| Selected ORF → BLAST | **Partly** | `build_orf_query` retains global coordinates; the remote path is untested (no live call) and the programmatic local path is blocked by D6. |
| Restriction → cloning | **Broken for circular DNA and ambiguity** | D3, D4, D5. The two modules have independent site finders with different semantics. |
| Primer ← restriction site selection | **Not wired** | `core/api.design_primers` takes only a target window; no restriction-site-driven primer path exists, although §10.2 of the Phase 6 brief describes one and the primer module has `restriction` references in its scoring. Needs a decision, not a fix. |
| VCF → parse → normalise → filter → statistics → export | **Works, with a silent failure mode** | Multi-allelic splitting, `Number=A/R/G` remapping, FILTER semantics and round-trip all pass their tests; D12 is the gap. |
| Mutation → coding consequences | **Works** | Start-loss restricted to codon 1, frame warnings, HGVS-*like* labelling all confirmed by the 44-test self-suite including the TEST_8 family. |
| Analysis → pipeline store → report | **Works** | JSON and Markdown exports produced; input SHA-256, parameters, warnings, limitations and per-analysis payload hashes present. |
| Literature | **Works offline** | Self-test passes on fixtures; live behaviour unverified (S11). |

---

## 8. Test coverage: what the 109 tests do not cover

| Test file | Tests |
|---|---|
| `test_p0_regressions.py` | 22 |
| `test_stage3_additions.py` | 21 |
| `test_stage1_additions.py` | 15 |
| `test_stage2_foundation.py` | 14 |
| `test_stage4_consolidation.py` | 10 |
| `test_primer_thermodynamics.py` | 7 |
| `test_stage5_protein_positions.py` | 5 |
| `test_main_smoke.py` | 2 |
| `test_builtin_selftests.py` | 1 (runs all twelve module self-tests) |
| `test_golden_master.py` | 12 golden scenarios |

Module references across all test files: orf 19, primer 14, mutation 11,
blast 11, variant 10, protein 8, restriction 7, cloning 7, literature 6,
**crispr 4, sequence 2, rna 2**.

Identified gaps:

1. **No test runs local BLAST**, although BLAST+ and a database are present.
2. **CRISPR is the least-tested large module** — 4 test references for 2,574
   lines; its 28 self-tests are its only real coverage, and they are the
   module's own fixtures, so D8–D10 pass unnoticed.
3. **No interface tests at all** — nothing checks menu structure,
   navigation keys, formatting or message wording (§5 of the brief requires
   these).
4. **No external-comparison tests** — every comparison in §3 was run as a
   throwaway script. None is in `tests/`, so none will catch a regression.
   These are exactly the tests that must skip gracefully when Biopython and
   primer3-py are absent.
5. **No test for the `cloning`-vs-`restriction` divergence** (D3–D5), the
   two `normalize_dna` contracts (D13), or alphabet detection (D14).
6. **No Python 3.10 run** and **no `pyright` run**, despite both being
   configured.
7. `tests/golden/orf_legacy_baseline.json` and the twelve snapshots protect
   behaviour but, by construction, cannot detect that the behaviour is wrong.

---

## 9. Documentation accuracy

`docs/AUDIT_STATUS.md`, `docs/CHANGES.md` and `docs/METHODS.md` are
unusually honest documents — they distinguish implemented, tested and
independently validated, and they name the validations that were not
performed. The audit found no overstated claim except:

| Document | Claim | Audit finding |
|---|---|---|
| `AUDIT_STATUS.md` P0-1 | Cites `blast_results/TEST_8_REFERENCE_blast.json` as evidence that coverage recomputes to 32.32 % | The *test* verifies this; the *file* still contains the defective 100 % values and an obsolete schema (D7). |
| `AUDIT_STATUS.md` P1-3 / P1-4 / P0-4 | "could not be run", "not possible in this environment" | All three are now closed by §3.1, §3.2, §3.4. The environment claim was wrong: the declared optional dependencies install cleanly. |
| `AUDIT_STATUS.md` P1-7 | Local BLAST "only testable with a database (intentionally excluded here)" | A test database is bundled in `data/local_test/` and BLAST+ is installed (§6). |
| `AUDIT_STATUS.md` P2-5 | ORF "Confidence" label "is now 'Heuristic score'" | True for the display only; the data model and all exports still say `confidence` / `High` / `Moderate` (D18). |
| `METHODS.md` | "Not implemented: 5′-partial ORFs without a start codon" | Correct, and worth keeping prominent: a 5′-truncated CDS is invisible to the default search. |
| `README.md` | "no third-party packages" | Accurate for the runtime. Should state that `[reference]` extras are needed for the external-comparison tests once those exist. |

`docs/METHODS.md` has no entry for the CRISPR scoring constants (S4) or for
the absence of a cut-site model (D8).

---

## 10. Prioritised implementation roadmap

Ordering principle: **scientific correctness before architecture before
presentation**, and within each tier, defects that silently produce wrong
numbers before defects that merely block a feature.

### Gate 0 — prerequisite, blocking all implementation — completed 2026-10-30
- **Git for Windows installed and verified** (2.56.0.windows.4). `git
  status`, `git diff` and `git commit` work; HEAD was `0dad3bb`; the two
  pre-existing working-tree differences were resolved by `5242fcb` and
  `5d8c73a` (§1). Gate 0 is a prerequisite only: Tier 1 still requires
  explicit approval before any application source file is changed.

### Tier 1 — silent scientific errors (do first)
1. **D12** variant contig matching and silent non-normalisation — a user is
   actively told that normalisation and REF checking happened when they did
   not. Add contig-name unification, and make an unmatched contig a
   warning on the result, not a counter nobody reads.
2. **D4 + D3 + D5** the cloning/restriction split — phantom sites from
   ambiguity, and circular DNA analysed as linear in the cloning planner.
   Converge both modules on one site finder with CONFIRMED/POSSIBLE status
   and a topology parameter.
3. **D2 + D1** ORF translation — make `translate_orf_sequence` read the
   ORF's own `genetic_code` and use the table's real start set. Then
   reconcile the ORF preview, the protein module and `core/api`.
4. **D14** alphabet detection — refuse, or at minimum flag, a non-DNA
   record at load time; set and check `SequenceRecord.alphabet`.
5. **D9** CRISPR on-target exclusion, and **D8** a cut-site model with a
   cited convention per system (or an explicit statement that none is
   modelled). D8 needs a scientific decision first, so it may follow.

### Tier 2 — misleading or dead functionality
6. **D7** regenerate or remove the stale BLAST artifact; if it is kept as a
   historical record, rename it and say so in the file.
7. **D6** fix the local-BLAST path gate (check for the database's component
   files, not the prefix), and **S7** populate `database_details` on that path.
8. **D10** either implement the five dead `DesignConfig` fields or delete
   them; **D11** default exports into the project's `exports/`.
9. **D15** decide the fate of `[11] AI Experiment Planner` — implement a
   first useful version, or label it explicitly as unavailable rather than
   "coming soon".
10. **D19** give `SynapsePipeline.run_orf` the full parameter set so a
    programmatic run can reproduce an interactive one.

### Tier 3 — tests that lock the gains in
11. Promote every §3 comparison into `tests/` as external-reference tests
    that **skip** when Biopython / primer3-py are missing, so the
    standard-library-only guarantee is preserved.
12. Add a local-BLAST integration test using `data/local_test/`, skipped
    when BLAST+ is absent.
13. Add regression tests for every Tier 1 and Tier 2 defect above.
14. Run `pyright` at `pythonVersion = 3.10` and install a real 3.10
    interpreter for a second suite run (**S12**, §1).
15. Add CRISPR coverage proportional to the module's size.

### Tier 4 — architecture
16. **D13** unify the two `normalize_dna` contracts; one validation module,
    one documented contract.
17. **D16** flatten the `_v5`–`_v11` primer ranking chain behind a
    behaviour-preserving test, as P2-1 recommends.
18. Extend `core/api.py` to the analyses it does not yet cover (cloning,
    variant filtering, mutation reporting, BLAST, literature) so the
    non-interactive layer is complete (**S14**).
19. Address **S9**: give the mutation module a hard, documented size limit
    with a clear refusal, instead of a note.

### Tier 5 — interface and language
20. Write `docs/TERMINAL_UI_STANDARD.md`, then implement shared header,
    section, menu, prompt, status and error renderers; retire the eleven
    separator widths and unify navigation keys (**D17**).
21. Write `docs/STYLE_GUIDE.md` and apply it, including **D18** — rename the
    `confidence` field and its `High`/`Moderate` classifications in the data
    model and exports, with a documented golden-master update.
22. Editorial pass over all user-facing text; add the missing METHODS
    entries for the CRISPR constants (**S4**) and the cut-site question (**D8**).
23. Fix the one language inconsistency found: `core/api.py:107` "analyse"
    against the project's "analyze".

### Still-open external validations
| Validation | Blocker | Route |
|---|---|---|
| `bcftools norm` comparison | No native Windows build | WSL or a Linux host; §3.7 is a hand-verified substitute |
| ORFfinder comparison | No local binary; web service | NCBI ORFfinder standalone, or a curated RefSeq CDS set |
| NEBcutter comparison | Web service | Superseded in substance by §3.2 (Biopython/REBASE on pUC19) |
| Live remote NCBI BLAST | Deliberately not run | One confirmed run, recording RID and provenance |
| Local BLAST on a real database | ~100 GB+ for `nt` | Not justified; `swissprot` (~0.5 GB) would be a proportionate middle step |
| Python 3.10 reproducibility | Interpreter absent | Install 3.10; run suite + pyright |

---

## 11. Audit limitations

1. **Git was unavailable while the audit was carried out.** History and
   refs were read from `.git` directly, and no diff-based review of the
   existing commits was possible while the findings were produced. Git
   2.56.0.windows.4 was installed afterwards, on 2026-10-30; §1 records
   what it verified, including the correction of the earlier "clean
   working tree" assumption.
2. **No application source file was modified**, per the standing
   instruction. All probes ran as throwaway scripts in the session
   scratchpad. This is now confirmed by Git: once it was available, the
   only tracked change in the working tree was the pre-existing
   `pyproject.toml` modification (§1). The Phase 6 files in the project are
   this document, `docs/FEATURE_STATUS.md`, and one added `.gitignore`
   line; the two housekeeping commits `5242fcb` and `5d8c73a` touch only
   `pyproject.toml` and `.gitignore`.
3. **Python 3.10 was not exercised** and **pyright was not run** (S12).
4. **No live network calls** were made to NCBI BLAST or to any literature
   provider. The pUC19 reference sequence was downloaded from NCBI
   E-utilities; nothing from the project was uploaded anywhere.
5. **Interactive menus were not driven end to end.** Every menu dispatcher
   was read and its entries inventoried, and the underlying functions were
   executed directly, but the terminal UIs were not exercised through
   simulated keystrokes beyond the existing `test_main_smoke.py`. Menu
   entries are therefore marked *reachable* in `FEATURE_STATUS.md`, not
   *execution-tested*, unless a test or a probe actually ran them.
6. **One audit task was interrupted.** While reading
   `analyses/crispr.py` lines 1199–1638 under Opus 5.5, the response was
   stopped by a safety classifier mid-turn. The file content was obtained,
   and the base-editing, prime-editing and ORF-integration sections of that
   range were subsequently audited under Opus 5 by targeted execution
   (D8–D10, S3, S5) rather than by re-reading the whole range. No finding
   was dropped, but that range was analysed by execution and targeted
   search rather than continuous reading.
7. **Skills and MCP providers:** the eleven project-local K-Dense skills
   under `.claude/skills/` were inventoried from disk. None was invoked;
   the external comparisons in §3 were run directly against Biopython and
   primer3-py, which is stronger evidence than a skill recommendation. The
   Bio Research MCP servers were not queried. `biorender`, `owkin`,
   `synapse` and `wiley` require authorisation and were unavailable.
8. **No claim in this document rests on an unexecuted assumption.** Every
   CONFIRMED finding quotes observed output; SUSPECTED and UNVERIFIED
   CONCERN items are labelled as such and state what would settle them.
