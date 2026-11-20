# Terminal interface standard

The shared visual and structural language of every Synapse terminal
interface. It is implemented in `ui/style.py` and enforced by
`tests/test_terminal_ui.py`, so it is a contract rather than advice.

Scope: the terminal interface only. Synapse has no graphical or web
interface, and Phase 6 does not add one.

The wider communication standard (terminology, units, scientific
writing, code documentation) is [`STYLE_GUIDE.md`](STYLE_GUIDE.md).

## 1. Why this exists

Before Stage F each component drew its own frame. A mechanical inventory
of `analyses/`, `ui/`, `core/` and `main.py` found:

| Problem | Measured |
|---|---|
| Separator widths in use | **11** (42, 45, 58, 60, 65, 68, 70, 72, 78, 88, 90) across 389 literal sites |
| Meanings of `[0]` | **7** (Back, Cancel, Exit, "Back / Exit", "Enter sequence manually", "Apply and return", "Choose another sequence") |
| Error prefixes | `ERROR:`, `Error:`, `Error ` |

BLAST framed its screens at 90 columns, literature at 88, CRISPR at 78,
most other modules at 70 with 45 for sub-blocks. The modules read as
separate programs (audit finding D17).

## 2. Layout

### Width

One width: `ui.style.RULE_WIDTH` = **70 columns**. It is the width already
most common in the code and it leaves margin inside an 80-column terminal.

Two rule characters:

- `=` (`major_rule()`) frames a screen.
- `-` (`minor_rule()`) divides a section inside a screen.

No module writes a separator literal. `major_rule()` and `minor_rule()`
are the only source.

**One deliberate exception.** A table rule must span its own table.
`analyses/orf.py` draws the candidate comparison with a rule that grows
with the number of columns and uses `RULE_WIDTH` as its floor. A rule that
is narrower than its table would be worse than an inconsistent one.

### Screen shape

```
======================================================================
MODULE OR SCREEN TITLE
======================================================================

SECTION TITLE
----------------------------------------------------------------------
Label:                 value
Label:                 value

GROUP
[1] First action
[2] Second action

[0] Back
----------------------------------------------------------------------
Selection:
```

Every part has one renderer in `ui.style`:

| Part | Renderer | Rule |
|---|---|---|
| Screen title | `header(title)` | Upper case, framed by two major rules, a blank line before. No "SYNAPSE" prefix: the title names the screen, and a module's main screen carries the name of its main-menu entry. |
| Context block | `field(label, value)` | Directly under the title. Labels padded to `FIELD_WIDTH` (22), so values start at column 24 in every module. A label is a short name; a sentence is a message, not a field. |
| Section | `section(title)` | Upper case over a minor rule, a blank line before. |
| Menu | `menu(groups, zero=...)` | The module's own groups and entries (5.3); each group after a blank line under an upper-case label; then a blank line, the `[0]` entry and a minor rule. |
| Prompt | `ask_choice()` / `ask(label, default=...)` | `Selection: ` for menus; a default in brackets. |
| Table row | `columns(values, widths)` | Values are converted to text before padding, so a column stays aligned whatever the number of digits. |
| End of a report | `print(major_rule())` | A report may close with a major rule. |

`field_line()` returns the same line as text, for reports that are built
as strings.

**Adoption.** Every interface file draws its screens through these
renderers (Stage F completion, 2026-11-02): `main.py`, `ui/menu.py`,
`ui/report.py` (which also draws the pipeline status that `core` hands it) and
the modules ORF, CRISPR, BLAST, literature, primer, variant analysis,
sequence overview, RNA, protein, restriction, cloning and mutation. What
they still print themselves is listed in section 9.

Since the closing review (2026-11-03) this also holds for the output the
rewriting passes had missed: validation and status messages go through the
message helpers, headings through `section()`, and option lists built in a
loop through `menu()`. A numbered result list that is not followed by a
choice is numbered `1.`, `2.` (literature results, BLAST evidence profile),
so `[n]` always marks a key the user can type.

**Tables** take their header and rows from one `columns()` definition;
values are left-aligned in their column and a table stays within the
70-column rule. **Text reports** that are built as line lists write
their key/value lines through `field_line()`, so a stored or exported
report aligns like the screen. Result rows that are lists of their own
(profile windows, restriction sites, CpG islands) keep their line form.

### Section order

The conceptual order is general information, input and context, main
analysis, advanced options, output, system. It is a guide, not a cage:
each module arranges the sections its science actually has. ORF analysis
organizes discovery, evaluation, comparison and selection; CRISPR design
organizes region, system, discovery, ranking, editing and specificity.
Forcing identical menus on different sciences would be worse than
inconsistency.

## 3. Navigation

**`[0]` always leaves the current screen.** That is the whole rule, and it
is what makes the interface safe to explore: a user can always press 0
without committing anything.

| Label | Where | Meaning |
|---|---|---|
| `Back` | any menu | Return to the menu above. |
| `Cancel` | a dialog that would otherwise commit an action | Abandon it; nothing is changed. |
| `Exit` | the main menu only | Leave Synapse. |

A screen that genuinely cannot be abandoned offers **no** `[0]` rather
than giving `[0]` a different meaning.

Letter keys name actions that are not navigation: `[R]` export report,
`[S]` select another sequence, `[P]` pipeline status, `[M]` enter a
sequence manually. Numeric keys are reserved for list items and ordered
options.

### The mutation sequence prompt

`analyses/mutation.py::_prompt_sequence` once offered no `[0]`: `[0]` had
meant "enter a sequence manually", which moved to `[M]`, and the four call
sites had no way to abandon the comparison. Since 2026-11-11 it offers
`[0] Cancel`, which raises `SequencePromptCancelled`; both menu branches
that ask for sequences catch it and return to the mutation menu without
comparing anything. The only screen without `[0]` is the sequence
selector, which cannot be left without choosing a sequence or `[F]`.

## 4. Input

| Helper | Use |
|---|---|
| `ask_choice()` | A menu selection. Renders `Selection: `. |
| `ask(label, default=...)` | A value. Pass `default=` only when it is the value the caller wants back for an empty answer. A hint such as `[none]` or `[50.0 mM]` stays in the label (`ask("From date [none]")`), and the caller maps the empty answer to its default. |
| `confirm(question, default_yes=False)` | A yes/no question, written `[y/N]` (or `[Y/n]` when yes is the default). `y`, `yes`, `j`, `ja` mean yes; anything else no, an empty answer the default. Every yes/no question in Synapse uses it. |
| `pause()` | `Press Enter to continue...` |

Invalid input never ends a workflow: the screen says `Invalid selection.`
and asks again.

## 5. Messages

One prefix per category, so a reader can tell at a glance what kind of
statement it is.

| Prefix | Meaning | Helper |
|---|---|---|
| `ERROR:` | The operation did not happen. | `style.error()` |
| `WARNING:` | It happened, but something needs attention. | `style.warning()` |
| `NOTE:` | Context worth knowing; nothing is wrong. | `style.note()` |
| `OK:` | Confirmation of a completed action. | `style.success()` |

| Situation | Category |
|---|---|
| An export or a save succeeded | `OK:` |
| An export, a file read or an input check failed | `ERROR:` |
| A view needs something first ("Run [1] BLAST Search first.", "No variants loaded.") | `NOTE:`, then `pause()` so it is read before the menu returns |
| An unknown menu key | `invalid_selection()`, then the same menu again; a dialog asks again rather than ending |

Scientific wording stays module-specific. A CRISPR off-target warning, a
primer thermodynamics warning and a variant normalization warning carry
the same prefix and entirely different text, because they say different
things. Uniform prefixes, individual science.

Note that guide quality flags such as `WARN: poly-T run` are **data**, not
interface messages. They are part of a result and keep their own
vocabulary.

## 6. Where the code lives

`ui/style.py` holds the presentation and nothing else: it decides how
something is drawn, never what it says or what it means.

`core/` draws nothing and imports nothing from `ui/` (Stage F,
recommendation 4, 2026-11-03). `PipelineContext.status_report()` returns
the pipeline status as data; `print_status(render)` hands that data to a
`StatusRenderer`, a callable the interface supplies. The terminal passes
`ui.report.print_pipeline_status`, which draws the screen with the shared
renderers. Until then `print_status()` drew the screen itself through
function-local imports of `ui.style`, which only deferred the dependency
to call time.

**The analysis modules do not import `ui` (P2-2, 2026-11-12).** The
dialogs of each module were moved, byte for byte, into
`ui/<module>_menu.py` (`ui/orf_menu.py`, `ui/crispr_menu.py`, ...): 130
definitions and 7 aliases, about 9,100 lines, from eleven modules. A
dialog is a definition that prints, reads input or uses a `ui.style`
renderer other than the string builders. The menu modules import the
science they call from `analyses.<module>`; `main.py` imports the dialogs
from `ui/`. Reports that the analysis modules build as strings
(`mutation.build_mutation_report`, `crispr.format_report`, the BLAST text
export, ...) take `field_line()`, `columns()` and the rules from
`core/textformat.py`, which `ui/style.py` re-exports, so a report and a
screen keep one layout. `SynapsePipeline.select_orf(chooser=...)` takes
the ORF dialog as an argument instead of importing it. An analysis module
imports `ui` only inside its `if __name__ == "__main__"` block, which
runs when the module is started as a script.

Evidence: the 55 rendered screens of `test_rendered_menus` and the first
screens of `test_terminal_ui` were recorded before and after the move and
are byte-identical; `LayeringTests` checks that no analysis module imports
`ui` outside a script block and runs every analysis of `core.api` in a
subprocess in which `ui` cannot be imported.

## 7. Enforcement

`tests/test_terminal_ui.py` checks the project, not a sample:

1. No separator literal outside `ui/style.py`, with the documented
   `orf.py` table exception.
2. Every `[0]` entry is `Back`, `Cancel` or `Exit`.
3. `Exit` appears only in the main menu and the sequence selector.
4. No `Error:`/`Warning:` prefix variants in user-facing messages.
5. `ui/style.py` imports nothing from `analyses/` or `core/`, so the
   presentation layer stays free of scientific logic. `core/` imports
   `ui` nowhere, not even inside a function; in an interpreter that
   refuses `ui`, every core module imports, the pipeline builds, reports
   and exports its status, and `python -m core.pipeline` runs
   (`LayeringTests`).
6. Menu numbers are contiguous, and every entry has a dispatch branch
   and every branch an entry (BLAST menu and main menu).
7. **Rendered screens** (`RenderedScreenTests`), driven with scripted
   input: all eleven menus (ten modules and the main menu) for frame,
   title (the main-menu name), context fields at the shared column,
   footer (blank line, `[0] Back` or `[0] Exit`, minor rule,
   `Selection: `) and the reply to an unknown key; the sequence overview
   (run through `main.py`) and the pipeline status for aligned fields and
   upper-case headings.
8. No interface file (all sixteen) prints a message prefix, menu entry,
   selection prompt or yes/no notation of its own, defines a rule or
   header helper, or titles a screen "SYNAPSE ...".
9. The variant menu follows the order of 5.3 (input, normalization,
   filtering, statistics, interpretation, reporting), with every entry
   dispatched.
10. **Direct output** (`DirectOutputTests`, read through the AST of all
    sixteen files): no option entry `[{i}] ...` printed in a loop over
    `enumerate()`; no validation or status sentence ("Invalid ...",
    "Please ...", "No ... supplied", "... cancelled.", "IMPORTANT:", ...)
    outside the message helpers; no upper-case heading, also not as the
    start of an f-string, outside `section()`; no heading drawn with its
    own dashes. The one named exception is the invalid-character list
    under the sequence validation error, which is a result line. On
    `21deb2a` the rule reports 68 sites. No message is drawn as a field
    (`style.field("Error", ...)`, `"... failed"`, `"... error"`); ten such
    lines in mutation, restriction and literature became ERROR and NOTE
    lines on 2026-11-11.
11. **Every menu rendered** (`tests/test_rendered_menus.py`): 40 dialogs
    driven with scripted keys; every `style.menu` call site of the
    sixteen files must be drawn by a scenario or a first-screen test
    (47 of 47); every drawn menu has unique keys, numbers 1..n, an
    allowed `[0]`, the entries it was given, a blank line before `[0]`
    and the closing rule; each dialog is checked for its entries, the
    reply to an unknown key, `0`, its error states and its result view.
    A missing file is reported as an ERROR by every dialog that reads one
    (sequence loader, VCF, CRISPR batch, BLAST query and database FASTA),
    and the dialog stays usable (`MissingFileTests`).

## 8. Numbering

Numbers follow the order of the work, not the order features were
added. When a new entry belongs in the middle of a group, the group is
renumbered rather than given the next free number at the end: the BLAST
menu's SEARCH block briefly read 1, 2, 3, 14 after the resume entry was
added, and now reads 1, 2, 3, 4.

Renumbering means rewriting the printed labels and the dispatch
branches from one mapping, and checking each branch against its new
label. In `analyses/blast.py` the rewrite must be scoped to
`blast_menu` by its AST bounds, because the file contains several other
numbered sub-dialogs that must not shift.

A change that breaks the standard fails the suite.

## 9. Open work

What remains outside the shared renderers, each for a stated reason:

- **Raw line input** (`input()` without a prompt) where a sequence or VCF
  is pasted line by line (CRISPR, variant analysis): the raw line must
  not be stripped or reformatted. `crispr._read_local_references` also
  reads raw lines but strips them; it belongs here because its prompt
  is the explanation printed above it, not a label.
- **Two CRISPR numeric helpers** take a complete prompt text from their
  callers (`input(prompt)`); converting them means changing every
  caller's prompt string.
- **Closing rules**: 18 direct `print(major_rule())` calls close a report
  or frame its closing summary line, which section 2 allows.
- **Result lines** of the analyses (profile windows, sites, islands,
  alignments, interpretation prose) remain module text: they are
  scientific content, not interface furniture. This includes result
  statements such as "No restriction sites detected." or "No ORFs
  found." (16 in ORF, CRISPR, protein, restriction and `ui/report.py`)
  and the invalid-character list under the sequence validation error.
- **Headings inside text reports** that are built as strings (the primer
  pair block's PRIMER PAIR #n, FORWARD, REVERSE, PAIR; the BLAST report's
  TOP HITS) are part of the stored report text, not printed screens.
- **The sequence selector** offers no `[0]` (section 3); it is drawn
  through `style.menu(zero=None)`. The cloning enzyme list does the same
  where a choice is required (`allow_back=False`).
- **Prompts of the converted option lists** keep their own label
  (`Database:`, `Program:`, `Record number:`), not `Selection:`; changing
  them was not part of recommendation 2.
- **Not yet covered by a rule**: the "RESULT:" lines of the module
  self-tests when a module is run as a script, and plain status lines such
  as "Loaded 2 FASTA record(s).".
  Resolved on 2026-11-17 (Phase 7): the statement below that no other
  hand-padded line remained was wrong; the single-primer analysis printed
  five ("ΔH:", "ΔS:", "ΔG at 25 °C:", "3' GC clamp:" twice). They use
  `style.field()` now, and so do three BLAST database-manager prompts that
  began with a carriage return. "Use .json, .csv or .txt." (BLAST export),
  "Cloning module could not be imported." and the literature date
  warning have their ERROR/WARNING prefix, "No hits were returned" is a
  NOTE. `tests/test_message_rules.py` fails on either form (the CRISPR
  dialog, outside Phase 7, is not checked).
  Resolved on 2026-11-13: the last two hand-padded label lines
  ("Retrieved before deduplication", the RNA "Input source") are drawn
  with `style.field()`; a search on 2026-11-13 found no other.
  Resolved on 2026-11-12: the literature limitation line is a NOTE; the
  BLAST hit selector accepts only the listed numbers and answers any other
  input with "Invalid selection." (it answered 2 of 1 with "ERROR: list
  index out of range" and showed a hit from the end of the list for -1);
  the BLAST query-type question answers an unknown key and asks again
  instead of cancelling.
- **Screens behind no menu test** (instrumented 2026-11-03): 46 of 101
  `header()` calls, 41 of 72 `ask()` calls and 2 of 7 `confirm()` calls
  are executed by the suite (line coverage of the dialogs since 2026-11-12:
  `docs/TEST_COVERAGE.md`). The rest are export paths, NCBI downloads
  and removals, remote submissions and similar dialogs that need files,
  tools or network; all 47 menus and all 33 `ask_choice()` calls run.
- **BLAST query type**: the explanation of an undetermined type is the
  exception text printed as it is.
