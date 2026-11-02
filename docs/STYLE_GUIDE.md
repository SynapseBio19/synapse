# Synapse style guide

How Synapse communicates: its scientific terms, its numbers, its messages,
its terminal screens and its code. This guide is the reference for all
future development. It describes the implementation as it stands, and the
rules that can be checked mechanically are checked by the test suite
(section 8), so the guide and the software cannot drift apart unnoticed.

The terminal layout has its own document,
[`TERMINAL_UI_STANDARD.md`](TERMINAL_UI_STANDARD.md); section 6 of this
guide summarizes it and does not repeat it.

## 1. Language

- **English** for every interface text, report, message, docstring and
  document. Synapse has no multilingual configuration; German or mixed
  text in the interface is a defect.
- **American spelling**, because the identifiers already use it
  (`normalize_sequence`, `analyze_orfs`): behavior, analyze, labeled,
  modeled, normalization, color. Words that are spelled with *-ise* in
  American English too (otherwise, precise, advertise, excise) are
  unaffected.
- **Sentence case** for prose and messages; screen and section titles are
  upper case, as the interface standard sets out.
- Write for a reader trained in molecular biology: precise, complete and
  plain. Academic register means accurate terms and stated limits, not
  long words.

## 2. Scientific terminology

### 2.1 Preferred terms

| Use | Not | Why |
|---|---|---|
| heuristic score (`heuristic_score`) | confidence, confidence score | The ORF score is a sum of fixed points. "Confidence" implies a probability. |
| evidence level (`evidence`: CONFIRMED, PROBABLE, AMBIGUOUS, INVALID) | confidence | Categories of evidence about the sequence alphabet, not probabilities. |
| heuristic evidence (`evidence`: low, moderate) | confidence | The signal-peptide screen is a sequence heuristic. |
| sequence type (DNA, RNA, PROTEIN) | alphabet (in prose) | The code keeps `alphabet` as a field name for compatibility; prose says what is meant. |
| ORF candidate | gene, protein | An ORF is a sequence feature. Coding potential is not shown by finding one. |
| modeled cut site | cut site | Synapse models where a nuclease cuts; it does not observe it. |
| predicted Tm | Tm, annealing temperature | A nearest-neighbor model value. It is not an annealing temperature. |
| sequence similarity (BLAST) | homology, identity of function | Homology and function are inferences that a hit alone does not establish. |
| off-target search in the supplied sequences | off-target analysis | The search never covers a genome unless one is supplied. |

**Confidence is reserved** for a statistical confidence (a confidence
interval or level). Synapse computes none, so the word appears only
in negations such as "not a statistical confidence". The test suite
rejects it as a data key.

### 2.2 Abbreviations

Defined once per screen or document where a reader may need it, then used
as is: ORF, CDS, UTR, PAM, sgRNA, PBS, RTT (prime editing), Tm, GC, NN
(nearest-neighbor), SNV, indel, VCF, HGVS, RID (BLAST request ID), HSP,
E-value, IUPAC.

### 2.3 Sequence and gene notation

- Sequences are upper case, 5'->3' unless a strand is named. Nucleotide
  ambiguity uses the IUPAC codes; amino acids the IUPAC-IUB one-letter
  codes, `*` for a stop.
- Strands are `+` and `-`. Reading frames are +1 to +3 and -1 to -3.
- Genetic codes are named by their NCBI table number: "NCBI table 11".
- Variants: HGVS-like protein changes with one-letter codes (`p.R12C`),
  not validated against the HGVS specification. HGVS writes a predicted,
  unconfirmed consequence in parentheses (`p.(R12C)`); Synapse does not do
  this yet, which is recorded as open work rather than claimed.
- Gene and protein symbols are reproduced exactly as given by the user or
  the source database. Synapse does not rename them, and a terminal cannot
  italicize them.

## 3. Numbers, units and coordinates

### 3.1 Units

One space between a number and its unit, except the percent sign.

| Quantity | Notation | Example |
|---|---|---|
| Percentage | number directly followed by `%` | `52.38%` |
| Temperature | `°C` after a space | `60.12 °C`, `ΔG at 25 °C` |
| Enthalpy, free energy | `kcal/mol` | `-128.80 kcal/mol` |
| Entropy | `cal/(mol·K)` | `-355.02 cal/(mol·K)` |
| Molecular mass | `Da` (`kDa` in the DNA sequence overview) | `12,345 Da` |
| Length | `nt` (single strand, primers, RNA), `bp` (duplex DNA, amplicons), `aa` (protein) | `20 nt`, `450 bp`, `149 aa` |
| Concentration | `mM` for ions, `nM` for oligonucleotides | `1.5 mM Mg²⁺`, `250 nM` |

Mass is shown in Da for proteins and nucleic acids alike. Da and g/mol are
numerically identical here, and one unit makes values comparable across
modules.

### 3.2 Precision

| Value | Format | Note |
|---|---|---|
| GC and other percentages | 2 decimals | |
| Tm, ΔH, ΔS, ΔG | 2 decimals | |
| BLAST identity, coverage, bit score | 2 decimals | |
| E-value | 3 significant figures (`1.23e-50`), `n/a` when absent | a configured threshold is shown exactly as set |
| Molecular mass | whole Da with thousands separator | kDa with 2 decimals |
| Heuristic scores | integer out of 100, or 2 decimals | always labeled heuristic |

Rounding is for display only. Results and exports keep full precision,
and no displayed value may be rounded so far that it changes an
interpretation (for example a ΔTm compared with a 2 °C limit).

### 3.3 Coordinates

- **Users see 1-based, inclusive positions**, in every menu and printed
  report. A range is written `117-135` and labeled `(1-based)` where the
  convention could be in doubt.
- **The code and exports keep their internal conventions**, which are
  documented per structure: `core.models.Interval` and primer `start`/`end`
  are 0-based half-open; ORF dicts are 0-based inclusive on the forward
  strand; VCF positions are 1-based. A field whose convention differs from
  its neighbors carries it in its name (`position_0based`,
  `site_position_1based`).
- Input from a user is converted at the menu, never inside the analysis.

## 4. Scientific writing

### 4.1 What kind of result is this?

Every reported value belongs to one of four kinds. The wording must make
the kind recognizable, and must never present a weaker kind as a stronger
one.

| Kind | Meaning | Synapse examples | Wording |
|---|---|---|---|
| **Measured** | Observed in an experiment | none: Synapse makes no measurements. Data read from files (sequences, VCF records, BLAST hits) are inputs, not Synapse measurements. | "the file reports", "NCBI returned" |
| **Calculated** | Follows deterministically from the input | GC content, molecular mass from composition, codon counts, ORF coordinates, restriction cut positions | plain statement: "GC content: 52.38%" |
| **Predicted** | Output of a model with assumptions | Tm (nearest-neighbor), hairpin and dimer screens, CRISPR on-target and modeled cut site, signal-peptide screen, Kozak context | "predicted", "estimated", "screen", with the model named |
| **Inferred** | A conclusion drawn from calculated or predicted values | coding consequence of a variant (assumes the CDS frame), homology from BLAST similarity, ORF ranking | "suggests", "is consistent with", "assuming ...", with the assumption named |

### 4.2 Structure of an interpretation

Where a module explains a result, it answers, in the order that fits the
science: what was analyzed; with which parameters; what was found; what it
means; what limits apply; whether experimental validation is needed. A
CRISPR off-target warning, a primer thermodynamics warning and a variant
normalization warning keep their own scientific content; only the
structure and the message prefix are shared.

### 4.3 Uncertainty and limits

- State limits as facts about the method: "Off-target search covers only
  the supplied sequence, not a genome."
- Negate over-claims explicitly where readers are likely to make them:
  "ORF detection alone does not prove biological expression."
- Conventions are called conventions. The ORF defaults, for example, are
  documented as NCBI ORFfinder screening conventions, not biological
  thresholds.
- Validation claims use the vocabulary of `docs/METHODS.md` and
  `docs/FEATURE_STATUS.md` and nothing stronger: known-answer test,
  implementation cross-check, compared with a named reference tool,
  execution-tested, not independently validated. "Validated" alone is not
  used; no Synapse result has been validated experimentally.

## 5. User communication

### 5.1 Message categories

One prefix per category, written through `ui/style.py`:

| Prefix | Meaning | Helper |
|---|---|---|
| `ERROR:` | The operation did not happen. | `style.error()` |
| `WARNING:` | It happened, but something needs attention. | `style.warning()` |
| `NOTE:` | Context worth knowing; nothing is wrong. | `style.note()` |
| `OK:` | A completed action, such as an export. | `style.success()` |

Guide and primer quality flags (`WARN: poly-T run`, `MARGINAL`) are result
data and keep their own vocabulary.

### 5.2 Recurring situations

These patterns are the standard. A message that departs from one is a
defect to correct, not an alternative style.

| Situation | Wording pattern |
|---|---|
| Unknown menu key | `Invalid selection.` (`style.invalid_selection()`), then the menu again. A range may follow: `Invalid selection. Please choose 0-6.` |
| Invalid value | Say what is required: `Positions must satisfy 1 <= start <= end <= 400.` |
| Missing input | `No DNA sequence loaded.` / `No DNA sequence selected.`: name the molecule the analysis needs. |
| Unsupported input | Name the type and the reason: restriction enzymes cut DNA, so RNA is refused. |
| Missing data in a result | `n/a` in tables, never `None` or an empty cell. |
| External service failure | Name the service and what happened, then what the user can do: `HTTP 429 from eutils.ncbi.nlm.nih.gov after 3 attempt(s); the source may be rate-limiting requests or temporarily unavailable.` |
| Missing dependency | Name the executable or database and where Synapse looked: `Local BLAST database 'x' is not installed or incomplete.` |
| Export, storage | `OK:` with the path written. |

Messages are statements, not exclamations, and contain no apology or
filler.

## 6. Terminal interface

Summarized from `TERMINAL_UI_STANDARD.md`, which is normative:

- One width, 70 columns: `=` frames a screen, `-` divides a section.
  Separators come only from `ui.style`.
- `[0]` always leaves the current screen: Back, Cancel, or Exit (main menu
  only). Letter keys name actions (`[R]` report, `[S]` select sequence).
- Menu numbers are contiguous and follow the order of the work.
- Key/value blocks align labels to 22 columns (`style.field`).
- Prompts: `Selection: ` for menus; defaults in brackets, used on empty
  input.

## 7. Code documentation

- **Python 3.10** compatibility and the standard library first; Biopython
  and primer3-py are test-only references.
- **Names**: `snake_case` functions and variables, `PascalCase` classes,
  `UPPER_CASE` constants. Private helpers start with `_`. Data keys use the
  terms of section 2 (`heuristic_score`, never `confidence`).
- **Docstrings**: one summary line in the imperative or descriptive mood,
  then what the reader cannot see from the signature: units, coordinate
  convention, assumptions, and what is not modeled. Code names in double
  backticks.
- **Type annotations** on public functions and dataclasses. New code adds
  no pyright diagnostics; the project baseline (134 at the start of
  Phase 6, 112 now) may only fall.
- **Comments** say why, not what: the reason for a convention, the source
  of a constant (with a PMID or a document), the audit item that
  motivated a change (`(audit D8)`).
- Scientific constants carry their source next to them, or are labeled as
  Synapse conventions.

## 8. Enforcement

The rules that can be checked are checked on every test run. A change
that breaks one fails the suite.

| Rule | Test |
|---|---|
| Separators, `[0]` meanings, message prefixes, contiguous menus, presentation layer free of science (section 6) | `tests/test_terminal_ui.py` |
| Rendered screens of the converted modules: frame, title, fields, footer, invalid-key reply; no bypass of `ui.style` (section 6) | `RenderedScreenTests` |
| "Confidence" is not used as a data key anywhere in the source (2.1) | `HeuristicScoreVocabularyTests` in `tests/test_phase6_tier2.py` |
| Unit notation in printed text: `%`, `°C`, `cal/(mol·K)`, `Da` (3.1) | `OutputNotationTests.test_printed_units_use_one_notation` |
| E-values with three significant figures (3.2) | `OutputNotationTests.test_blast_text_report_shows_e_values_with_three_significant_figures` |
| Primer positions shown and read 1-based (3.3) | `OutputNotationTests` (two tests) |
| American spelling in code and documentation (1) | `SpellingStandardTests` |
| Python 3.10 compatibility: 3.10 grammar, no f-string that needs 3.12 (7) | `Python310CompatibilityTests` |
| A result's own warnings reach the exported summary (4.3, 5) | `ExportedWarningsTests` |

Rules that need judgment (sections 2.3, 4 and 5.2) are applied in review.
When a new module or message is added, it follows this guide; when the
guide is wrong, the guide is changed in the same commit as the code.
