# Stage F re-examination: the terminal interface as implemented

Requested by the owner on 2026-11-02. The question: was Stage F completed,
or were mainly separators, prefixes and some navigation conventions
standardized? This review changes nothing. Every figure below was
measured on `fc13aad`:

- source scans of `main.py`, `analyses/`, `ui/` and `core/pipeline.py`
  (AST and text);
- rendered screens: the first screen of each module menu, driven with
  scripted input in a scratch process.

## Verdict

**Stage F delivered a common formatting basis, not a unified user
interface.**

What was standardized and is enforced:

- one rule width (70) and two rule characters through
  `ui.style.major_rule()`/`minor_rule()`;
- the meanings of `[0]` (Back, Cancel, Exit);
- the spelling of the `ERROR:`/`WARNING:`/`NOTE:` prefixes where a
  prefix is used;
- the numbering of the BLAST and main menus;
- the text "Invalid selection." in most menus.

What was not done is the larger part of Master Prompt sections 5.1, 5.4
and 5.5. `ui/style.py` contains shared renderers for headers, sections,
fields, menu options, prompts, confirmation, pauses and every message
category, but **no module calls any of them**. Every module still
renders its screens with its own `print` and `input` calls (about 2,500
`print` and 140 `input` calls in the files reviewed). Headings,
spacing, alignment, prompts, context blocks, result tables, message
categories and precondition handling still differ between modules,
including the six named by the owner.

The Stage F documents describe the helpers as the standard
(`TERMINAL_UI_STANDARD.md`: "`ui.style` provides `header()`, `section()`
and `field()`"; "`ask_choice()` ... Renders `Selection: `"). The final
report says "one visual language across all terminal screens" and, for
RNA and protein, "shared frames, prompts and messages". These statements
overstate the implementation.

---

## 1. Which modules use `ui/style.py`

Calls of `ui.style` functions per file (AST count; `dataclasses.field`
excluded):

| File | `print` | `input` | `major_rule` | `minor_rule` | any other `ui.style` helper |
|---|---:|---:|---:|---:|---|
| `main.py` | 90 | 14 | 17 | 0 | none |
| `analyses/orf.py` | 233 | 9 | 6 | 22 | none |
| `analyses/crispr.py` | 189 | 19 | 18 | 0 (own `_line()`) | none |
| `analyses/blast.py` | 269 | 30 | 2 | 3 (own `_line()`) | none |
| `analyses/literature.py` | 83 | 11 | 4 | 1 | none |
| `analyses/primer.py` | 185 | 10 | 21 | 11 | none |
| `analyses/variant_analysis.py` | 135 | 16 | 17 | 14 | none |
| `analyses/mutation.py` | 107 | 10 | 18 | 15 | none |
| `analyses/cloning.py` | 193 | 12 | 15 | 18 | none |
| `analyses/restriction.py` | 140 | 5 | 17 | 18 | none |
| `analyses/rna.py` | 287 | 2 | 19 | 37 | none |
| `analyses/protein.py` | 372 | 2 | 8 | 49 | none |
| `analyses/sequence.py` | 7 | 0 | 2 | 1 | none |
| `ui/menu.py` | 30 | 2 | 4 | 4 | none |
| `ui/report.py` | 158 | 0 | 3 | 19 | none |
| `core/pipeline.py` | 28 | 0 | 5 | 2 | none |

`header`, `section`, `field`, `option`, `options`, `back`, `cancel`,
`exit_option`, `ask`, `ask_choice`, `pause`, `confirm`, `error`,
`warning`, `note`, `success`, `invalid_selection` and `notes` have **zero
callers** outside `ui/style.py`. The only exceptions are its own tests
(`test_helpers_render_the_documented_prefixes`).

## 2. Module-local formatting and menus

Every module draws its own screens. Local helpers that duplicate or
bypass the shared layer:

| File | Local helper | Note |
|---|---|---|
| `analyses/blast.py:3807` | `_line(char, width=RULE_WIDTH)` | a second rule function; with the one in CRISPR, 57 call sites |
| `analyses/crispr.py:1929` | `_line(char, width=RULE_WIDTH)` | same |
| `analyses/orf.py:2936` | `_print_orf_header(title)` | a module-own header renderer |
| `analyses/blast.py:4690` | `print("-----------------------------------------------")` | a literal 47-character rule; the Stage F test does not detect literal runs |
| `analyses/blast.py:3813, 3815, 3845` | `print("\r\n" + _line("="))` | prints an explicit carriage return, unlike every other module |

Module-local function names that render or prompt (pattern `print`,
`menu`, `_ui`, `display`, `show`, `header`, `format`, `table`): cloning
19, crispr 13, protein 13, orf 12, primer 12, restriction 10, rna 10,
variant 10, blast 7, mutation 5, literature 3.

## 3. Headings, spacing, alignment, options and messages

### Screen titles

Collected from the text between two `major_rule()` lines:

- **Case.** Upper case almost everywhere, but the main menu's title is
  "Choose Analysis" (`ui/menu.py:6`).
- **Prefix.** Three forms in use:
  - `SYNAPSE — ` with an em dash (5 titles), e.g. `analyses/primer.py:897`
    "SYNAPSE — PRIMER DESIGN" and `analyses/restriction.py:2131`;
  - `SYNAPSE - ` with a hyphen (2 titles), e.g. `main.py:127`
    "SYNAPSE - LOAD SEQUENCE FILE";
  - no prefix (46 titles).
- **Names differ between the main menu and the module.**

| Main menu entry | Module title |
|---|---|
| `[8] CRISPR Guide Design` | "CRISPR ANALYSIS ENGINE" (`analyses/crispr.py:2564`) |
| `[10] Literature Search` | "Synapse Literature Intelligence v2.1.1" (`analyses/literature.py:2191`) |
| `[6] Primer Design` | "PRIMER DESIGN" in the menu, "SYNAPSE — PRIMER DESIGN" on the result screen |

### Section headings

Taken from the line printed before a `minor_rule()`:

| Module | Style |
|---|---|
| ORF | Title Case ("ORF Summary", "Reading Frame Statistics") |
| RNA, restriction, `ui/report.py` | Title Case |
| Cloning | UPPER CASE ("DESIGN EVALUATION", "VECTOR") |
| Primer | mostly UPPER ("TARGET REGION", "TM MODEL"), 3 Title Case |
| Protein | **both within one module**: 26 UPPER, 20 Title Case |
| CRISPR, BLAST, variant menus | group labels in UPPER ("GUIDE DESIGN", "RESULT ANALYSIS") without a rule |

### Key/value alignment

`ui.style.FIELD_WIDTH` is 22, and no module uses it. Measured as the
length of "Label:" plus padding in f-strings, most frequent widths
first:

| Module | Value column |
|---|---|
| protein, RNA, `ui/report.py` | 21 |
| ORF | 23, 21, 25 |
| primer | 23, 24, 28 |
| CRISPR | 26, 23, 17 (20 distinct widths) |
| BLAST | 18, 7, 20 (21 distinct widths) |
| variant | 19, 20 |
| mutation | 18, 24, 19 |
| literature | 9, 5, 23 (17 distinct widths) |

The context blocks of the menus show it directly:

- ORF: `Sequence length: 900 bp`
- CRISPR: `Sequence: demo` / `Length:   900 bp`
- BLAST: `Query context:    demo`
- primer: `Input:  Sequence` / `Length: 900 bp`
- variant: a blank line, then `Loaded variants: 0`

### Menu prompts

| Prompt | Uses | Where |
|---|---:|---|
| `Selection: ` | 24 | blast, cloning, crispr, literature, menu, mutation, primer, protein, restriction, rna, variant |
| `Select option: ` | 2 | `analyses/orf.py:3247` (settings) and `:3375` (**the ORF main menu**) |
| `Selection [1]: ` | 2 | `analyses/primer.py:3882, 3931` |
| `Select an option: ` | 1 | `main.py:291` (protein record menu) |
| `Select ORF: ` | 1 | `analyses/orf.py:3180` |
| `Select a sequence: ` | 1 | `ui/menu.py:58` |
| `{…} selection: ` | 1 | `analyses/mutation.py` |

Confirmation prompts use `[y/N]` (5) and `(y/N)` (1,
`analyses/literature.py:2291`); `style.confirm()` is unused.

### Menu layout

Rendered first screens:

| Module | Framed title | Context block | Option groups | Blank line before `[0]` | Rule before the prompt | Prompt |
|---|---|---|---|---|---|---|
| ORF | yes | yes | none | no | **no** | **`Select option:`** |
| CRISPR | yes | yes | 6 labeled groups | no | yes | `Selection:` |
| BLAST | yes | yes | 4 labeled groups | no | yes | `Selection:` |
| Literature | **no frame**; banner "Synapse Literature Intelligence v2.1.1" | no | 1 label | no | **no** | `Selection:` |
| Primer | yes | yes | none | yes | yes | `Selection:` |
| Variant | yes | blank line, then count | 4 labeled groups | yes | yes | `Selection:` |
| Main menu | yes ("Choose Analysis") | no | groups separated by rules, not labels | yes | no (blank line) | `Selection:` |

Different groupings are legitimate (5.3). The frame, context block,
spacing, footer and prompt around them are not: those are exactly what
5.1 and 5.2 ask to be uniform.

### Status messages

Literal prefix use per module:

| Module | `ERROR:` | `WARNING:` | `NOTE:` | `OK:` |
|---|---:|---:|---:|---:|
| blast | 14 | 0 | 2 | 0 |
| primer | 10 | 0 | 0 | 0 |
| variant | 8 | 0 | 0 | 0 |
| crispr | 7 | 1 | 0 | 0 |
| restriction | 0 | 1 | 0 | 0 |
| orf, literature, mutation, cloning, rna, protein, main | 0 | 0 | 0 | 0 |

The same categories appear without a prefix elsewhere:

- **Failures without `ERROR:`**:
  - `print(f"Export failed: {exc}")` in `analyses/blast.py:4747`,
    `analyses/crispr.py:2541, 2558` and `analyses/literature.py:2356`;
  - `print("Export failed: permission denied.")` in
    `analyses/crispr.py:2538, 2555`;
  - `print(f"Not available: {exc}")` in the ORF settings.
- **Successes without `OK:`**: `OK:` is used nowhere. Examples are
  `Exported: {path}` (`analyses/blast.py:4745`, `crispr.py:2536`,
  `literature.py:2352`), `VCF exported: {filename}`
  (`variant_analysis.py:2750`) and `Report written: ...`.
- **Precondition messages without a category**: "Run [1] BLAST Search
  first." (`analyses/blast.py:4612, 4644`); "No variants loaded." (9
  sites, e.g. `variant_analysis.py:2607`); "No ORF candidates
  available."; "No DNA sequence selected.".

The invalid-selection reply, rendered by entering `x`:

| Module | Reply |
|---|---|
| ORF | `Invalid selection. Please choose 0-6.` |
| CRISPR, BLAST, literature, variant | `Invalid selection.` |
| Primer | `ERROR: Invalid selection.` (the dialog ends) |
| Main menu | `Invalid selection.` (`main.py:979`, from the session loop; `ui/menu.show_menu` only returns the input) |

After a precondition message, BLAST and variant redraw the menu at once
(rendered), and CRISPR does the same in the code (`print(...)` followed
by `continue`). Result screens elsewhere end with "Press Enter to
continue..." (rendered for the ORF candidate list). Whether a message
stays readable before the menu returns therefore depends on the module
(5.5, "Input retry"). The pause behavior of every module was not
rendered for this review.

## 4. Result views and error presentation

Result screens are module-local throughout. One example from the ORF
candidate list, rendered:

```
ID        Frame    Length     Protein    Start   Status                  Conf.
----------------------------------------------------------------------
ORF-001  +3      207 bp     68 aa     ATG     COMPLETE                80/100
ORF-008  -3      87 bp     28 aa     ATG     COMPLETE                78/100
```

- **The header and the rows do not align.** The header pads `ID` to 10
  columns, the rows print `ORF-001` plus two spaces (9). Length and
  protein length are printed at their natural width plus a fixed suffix,
  so every following column shifts when the number of digits changes
  (`analyses/orf.py:2991-3008`).
- **The column is still labeled `Conf.`** (`analyses/orf.py:2997`) for
  the value D18 renamed to heuristic score. This is a D18 remnant in a
  display string that the D18 checks did not cover: they look at data
  keys, not column labels.
- The status column prints `COMPLETE` here, while the summary above it
  prints `Status: Complete`.

The other result views (BLAST hit list, CRISPR guide table, primer pair
summary, variant table) each have their own column widths, headers and
number formats. Only the unit notation (Stage G) is shared and tested.

Error presentation follows no shared path (section 3, "Status messages").
`style.error()`, `style.warning()`, `style.note()` and `style.success()`
exist and are tested in isolation, but nothing calls them.

## 5. What the Stage F tests check

`tests/test_terminal_ui.py` (12 tests):

| Test | Checks | Does not check |
|---|---|---|
| `test_no_module_draws_its_own_separator` | no `"=" * n` / `"-" * n` outside `ui/style.py` | literal rule strings (`blast.py:4690`), local wrappers (`_line`) |
| `test_one_width_for_the_whole_application` | `RULE_WIDTH == 70` and rule lengths | - |
| `test_every_ui_file_uses_the_shared_rules` | at least 15 files contain `major_rule()` or `minor_rule()` | whether a screen is framed at all (literature) |
| `test_zero_always_leaves_the_screen` | every `print("[0] ...")` label is Back, Cancel or Exit | blank-line placement, `[0]` rendering through `style.back()` |
| `test_exit_belongs_to_the_main_menu_only` | where `Exit` appears | - |
| `test_is_back_recognises_the_key` | `style.is_back` | that menus use it |
| `test_blast_menu_is_contiguous_and_fully_dispatched`, `test_the_main_menu_is_contiguous` | numbering in two menus | the other menus |
| `test_message_prefixes_are_uniform` | no `print("Error ...` / `Warning:` / `Note:` variants | whether failures and successes carry a prefix at all |
| `test_helpers_render_the_documented_prefixes` | the helpers' own output | that any module uses them |
| `test_style_module_carries_no_scientific_logic`, `test_core_has_no_module_level_dependency_on_ui` | layering | - |

**The tests check separators, `[0]` labels, prefix spelling, two menus'
numbering and layering.** Nothing tests:

- headings, title format or case;
- field alignment;
- prompts;
- spacing and footer;
- result-table alignment;
- the presence of message categories;
- precondition behavior;
- rendered output.

There is no snapshot of a rendered screen. The tests are honest about
what they check, but what they check is a small part of 5.1.

(Minor: the test name `test_is_back_recognises_the_key` uses British
spelling; the spelling test's word list does not include "recognises".)

## 6. The six named modules side by side

| Aspect | ORF | CRISPR | BLAST | Literature | Primer | Variant |
|---|---|---|---|---|---|---|
| Title | `ORF ANALYSIS` (own `_print_orf_header`) | `CRISPR ANALYSIS ENGINE` (name differs from the main menu) | `BLAST SEARCH` | unframed banner with version | `PRIMER DESIGN`; results `SYNAPSE — PRIMER DESIGN` | `VARIANT ANALYSIS` |
| Context block | length, count, search, best candidate (17-col labels) | sequence, length, region, system (10-col) | context, source, length (18-col) | none | input, length (8-col) | blank line + count |
| Grouping (5.3) | flat list; settings is entry 6 | 6 groups: INPUT / CONTEXT, GUIDE DESIGN, EDITING, SPECIFICITY, OUTPUT, SYSTEM | 4 groups: SEARCH, RESULT ANALYSIS, OUTPUT, SYSTEM | one label | flat, 2 entries | 4 groups: INPUT, ANALYSIS, OUTPUT, SYSTEM |
| Footer and prompt | no rule; `Select option:` | rule; `Selection:` | rule; `Selection:` | no rule; `Selection:` | blank, rule; `Selection:` | blank, rule; `Selection:` |
| Invalid key | `Invalid selection. Please choose 0-6.` | `Invalid selection.` | `Invalid selection.` | `Invalid selection.` | `ERROR: Invalid selection.`, dialog ends | `Invalid selection.` |
| Missing precondition | `No ORF candidates available.` (not rendered) | `No DNA sequence selected.`, menu redrawn (code) | `Run [1] BLAST Search first.`, menu redrawn (rendered) | - | - | `No variants loaded.`, menu redrawn (rendered) |
| Rule helper | `major_rule`/`minor_rule` | own `_line()` + `major_rule` | own `_line()`, one literal rule | `major_rule` in sub-screens only | shared | shared |
| Messages | no prefixes | `ERROR:` x7, `Export failed:` unprefixed | `ERROR:` x14, `Export failed:` unprefixed | none; `Export failed:` unprefixed | `ERROR:` x10 | `ERROR:` x8; exports unprefixed |

The groupings follow 5.3 reasonably for CRISPR, BLAST and variant. The
ORF and primer menus are flat; 5.3 asks for ORF settings and the primer
stages (design parameters, thermodynamic evaluation, specificity,
pairing, output) to be organized, and the primer menu exposes only two
entries.

## 7. Against the Master Prompt

| Requirement | State |
|---|---|
| 5.1 One specification covering header formatting, titles, section hierarchy and naming, sequence/region/strand information, context sections, result presentation, spacing, alignment | **Partly.** `TERMINAL_UI_STANDARD.md` covers width, rule characters, screen shape, `[0]`, prompts and prefixes. It sets no rules for section-heading case, title prefix, context-block content and order, result tables, footer and blank lines. Where it does set a rule (field width 22, `Selection:`), the code does not follow it everywhere. |
| 5.2 Shared structural template | **Partly.** Followed by CRISPR, BLAST and variant (groups, OUTPUT/SYSTEM); not by ORF, primer, literature or the main menu. |
| 5.3 Scientific individuality | **Met** where menus are grouped; ORF and primer organization is thinner than 5.3 describes. |
| 5.4 Shared header, section, menu, prompt, navigation, status, warning, error, result and export renderers | **Not met.** The renderers exist but are unused, so "a shared change applicable consistently across modules" is only true for rule width. |
| 5.5 Navigation: back, cancel, return, invalid selection, exit, confirmation, input retry | **Partly.** `[0]` meanings are uniform. The invalid-selection reply differs (ORF with a range, primer with `ERROR:` and an ended dialog); retry behavior differs (pause vs redraw); confirmation has two notations; there is no "return to main menu". |
| 5.6 Complete interface audit and `TERMINAL_UI_STANDARD.md` | **Partly.** The document exists, but the audit measured separators, `[0]` meanings and prefixes only; it records three deviations, while this review finds the larger set above. |
| 6.1 Consistent candidate tables | **Not met** for ORF (misaligned, `Conf.`); other tables are module-specific. |
| 16 "All interfaces follow one shared design language" | **Not met** beyond rules and `[0]`. |

## 8. Remaining inconsistencies

1. No module uses the shared renderers (`header`, `section`, `field`,
   `option(s)`, `back`, `ask_choice`, `confirm`, `pause`,
   `error`/`warning`/`note`/`success`, `invalid_selection`).
2. Local rule helpers in BLAST and CRISPR (`_line`, 57 calls), one literal
   rule (`blast.py:4690`), and explicit `"\r\n"` in BLAST headers.
3. Screen titles: three prefix forms (em dash, hyphen, none), one Title
   Case title, and module titles that differ from the main-menu entries
   (CRISPR, literature).
4. Section headings: UPPER and Title Case across modules and inside
   `protein.py`.
5. Key/value value columns from 18 to 28 characters; `FIELD_WIDTH` unused.
6. Prompts: `Select option:` (ORF main menu), `Select an option:`
   (protein record menu), `Selection [1]:` (primer), alongside
   `Selection:`.
7. Menu frame: literature is unframed; footer rule and the blank line
   before `[0]` vary per module.
8. Messages: `OK:` is never used; failures and precondition notes appear
   without a category; exports are reported in four wordings.
9. Invalid selection: ranged in ORF, `ERROR:` with an ended dialog in
   primer; redraw vs pause after precondition messages.
10. Confirmation: `[y/N]` vs `(y/N)`.
11. ORF candidate table: misaligned columns, the `Conf.` label (a D18
    remnant), and `COMPLETE` vs `Complete`.
12. Result tables in BLAST, CRISPR, primer and variant follow no shared
    table convention.
13. Tests cover none of points 1 and 3-12, and detect point 2 only in
    the `"-" * n` form.

## 9. What a completion of Stage F would involve

For the owner's decision. Nothing here was started.

1. **Adopt the existing helpers** module by module, beginning with the
   six named ones: `style.header()` for every screen title, `section()`
   for headings, `field()` for context blocks, `options()`/`back()` and
   `ask_choice()` for menus, and the message helpers for every failure,
   precondition and success. This is a presentation change: the golden
   master (data only) is unaffected. The end-to-end menu tests that
   match printed text would need review.
2. **Extend the standard** where it is silent: title format (one prefix
   rule), heading case, context-block content and order (sequence,
   length, region, strand, context, as 5.2 lists), footer and spacing,
   table conventions, the precondition message pattern, and retry
   behavior.
3. **Remove the local rule helpers and the literal rule.**
4. **Fix the ORF table**: aligned columns, and the `Conf.` label renamed
   to the heuristic score.
5. **One invalid-selection behavior**: the same reply and the same retry
   in every menu, through `style.invalid_selection()`.
6. **Test the rendered interface**, not only source patterns: render
   each module's first screen with scripted input and assert the frame,
   title, footer, prompt and alignment; and fail on any `print` of a
   heading or message that bypasses `ui.style`.
7. **Correct the documents** that describe the helpers as the implemented
   standard (`TERMINAL_UI_STANDARD.md` sections 2 and 4, the final
   report's "one visual language" and its section E rows).

---

## 10. Follow-up: the six modules converted (2026-11-02)

On the owner's instruction ("start with six"), commits `2e47ff9`, `e7eaf58`, `5cb4c9a`, `b181ae3`, `fb4070a`, `fdc968d`, `b80fd37`
route the ORF, CRISPR, BLAST, literature, primer and variant screens
through `ui/style.py`:

| Finding above | State in the six modules |
|---|---|
| 1. Shared renderers unused | `header`, `section`, `field`, `menu`, `ask_choice`, `ask`, `pause`, `error`, `warning`, `note`, `success` and `invalid_selection` are used throughout; a test fails if a converted module prints a prefix, menu entry or selection prompt of its own. |
| 2. Local rule helpers, literal rule, `"\r\n"` | Removed (`_line` in BLAST and CRISPR, `_print_orf_header`, `_primer_ui_header`). |
| 3. Title prefixes and names | No "SYNAPSE" prefix; CRISPR is "CRISPR GUIDE DESIGN" and literature "LITERATURE SEARCH", the names of their main-menu entries. |
| 4. Heading case | Sections upper case through `section()`. |
| 5. Field alignment | Values at column 24 in all six. |
| 6. Prompts | `Selection:` everywhere; defaults through `ask()`. |
| 7. Menu frame | Literature framed; one footer in all six. |
| 8. Messages | Prefixed by category; exports `OK:`/`ERROR:`; preconditions `NOTE:` with a pause. |
| 9. Invalid selection | One reply; the primer Tm dialog asks again. |
| 10. Confirmation | Unchanged (`[y/N]` and one `(y/N)` in literature). |
| 11. ORF candidate table | Aligned through `columns()`, labeled "Heuristic". |
| 12. Other result tables | Unchanged. |
| 13. Tests | `RenderedScreenTests` render the six first screens and check frame, title, fields, footer and the invalid-key reply; they fail on the code before the conversion. |

Still open: the other interface files (`TERMINAL_UI_STANDARD.md` section 9),
the variant menu order, other result tables, dialog value prompts and the
confirmation notation.

## 11. Follow-up: the remaining ten files converted (2026-11-02)

On the owner's instruction the other ten interface files were converted
(`1ac77a4`, `bf6ae87`, `9b7f2ef`, `85c22ac`, `a4df3d4`, `3853d87`, `69e3f0f`, `7970e5a`, `7ecb88d`): `main.py`, `ui/menu.py`, `ui/report.py`, sequence overview,
RNA, protein, restriction, cloning, mutation and `core/pipeline.py`.

| Finding | State in all sixteen interface files |
|---|---|
| 1. Shared renderers unused | All files use them: 779 `field`, 172 `section`, 142 `field_line`, 102 `header`, 36 `menu`, 74 `ask`, 31 `ask_choice`, 58 `pause`, 7 `confirm`, 17 `columns`, and 168 message helpers. |
| 2. Local helpers | Also removed: `_cloning_ui_print_header`/`_footer`, `_protein_ui_header`. |
| 3. Titles | Each menu carries its main-menu name ("MAIN MENU", "RNA ANALYSIS", "PROTEIN ANALYSIS", "RESTRICTION & CLONING ANALYSIS", ...); no "SYNAPSE" prefix anywhere. |
| 4./5. Headings, fields | Upper-case sections and values at column 24 in the converted screens (19 headings were still printed directly, see section 12), including the sequence overview, the pipeline status and the text reports of CRISPR, BLAST, primer, variant, literature and mutation. |
| 6. Prompts | `Selection:` in every menu; 52 value prompts through `ask()`. |
| 7. Menu frame | One footer everywhere; the sequence selector and the mutation sequence prompt, which cannot be left, offer no `[0]` through `menu(zero=None)`. |
| 8. Messages | Categories for the messages the passes recognized (31 uncategorized ones remained, see section 12); the checkmark lines ("✓ ...") are `OK:` messages. |
| 9. Invalid selection | One reply in every menu ("Please choose 0-5" and similar variants removed). |
| 10. Confirmation | All six yes/no questions through `confirm()` (`[y/N]`/`[Y/n]`). |
| 11./12. Tables | One `columns()` definition per table. Four tables were misaligned and are fixed: BLAST database status (Status 42 vs 24 columns), CRISPR guide list (position column), both protein composition tables (3L column); BLAST hits, restriction enzymes and taxonomy now fit the 70-column rule. |
| 13. Tests | `RenderedScreenTests` render the first screen of all eleven menus, the sequence overview and the pipeline status; the source rule covers all sixteen files; a numbering test checks the variant order. They give 11 failures on the code before this round. |

Variant menu (requirement 4 of the owner's instruction): reordered to the
5.3 sequence (normalization [3] and filtering [4] before statistics [5]
and interpretation [6]); every dispatch branch was checked against its
new label.

Also corrected on the way:

- **A defect from the first round.** Four BLAST and variant yes/no
  questions had been turned into `ask(..., default='y/N')`. The answer
  was still read correctly, but the code treated the notation as a value.
- **A message.** The ORF hand-over said the RNA codon view uses the
  standard code, which has been untrue since P1-9.

**Verdict as first written, corrected by the closing review.** This
paragraph originally said that Stage F was "implemented as a unified
interface" and that "every screen, menu, prompt, message, table and
text report" went through `ui/style.py`. That overstated the result.

The independent closing review (`docs/STAGE_F_FINAL_REVIEW.md`,
2026-11-03) measured the following:

- **Adoption.** All sixteen files use the shared renderers. But 59
  direct prints remained, and none of them was documented:
  - 31 validation or status messages without a category;
  - 19 upper-case headings printed directly;
  - 9 dynamic option lists (`print(f"[{i}] ...")`) in BLAST and CRISPR.
- **Test coverage.** The rendered tests drive the first screen of every
  menu, which is 17 of the 38 `style.menu` call sites. The source rule
  matched literal keys only.
- **Behavior.** It was not preserved from round 1 until `1516339`
  (section 12 below).
- **`core/pipeline.py`.** It still depends on `ui/` at call time. The
  function-local imports only defer the dependency.

The state after this round is a broad adoption of the shared renderers,
not yet a unified interface. Section 12 records what followed.

## 12. Enter regression and closing review (2026-11-03)

**Enter regression (`1516339`).** The rewriting pass had turned
`Label [hint]: ` prompts into `ask(label, default=hint)`, so Enter
returned the hint text instead of an empty answer. Fifteen prompts were
affected; four dialogs broke:

- the literature dates;
- the primer conditions;
- the BLAST word size and filter;
- the cloning strategy count, which looped endlessly.

The defect came in with `5cb4c9a`, `b181ae3`, `fb4070a` and `85c22ac`.
Until `1516339` the claim "the workflow is preserved" was false for
those dialogs. `EmptyAnswerTests` now drive them with an empty answer.

**Closing review.** `docs/STAGE_F_FINAL_REVIEW.md` checks the history
`623a0a8..1516339` and records the open points:

- the 59 direct prints;
- the deferred core→ui dependency;
- the stale menu numbers in `FEATURE_STATUS.md`;
- the test coverage (17 of 38 menus).

The owner released its recommendations 1–3: correct the claims,
convert the 59 prints, and correct `FEATURE_STATUS.md`. Recommendations
4 (core/ui) and 5 (full menu coverage) are deferred.
