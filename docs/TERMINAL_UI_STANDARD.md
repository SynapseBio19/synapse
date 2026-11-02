# Terminal interface standard

The shared visual and structural language of every Synapse terminal
interface. It is implemented in `ui/style.py` and enforced by
`tests/test_terminal_ui.py`, so it is a contract rather than advice.

Scope: the terminal interface only. Synapse has no graphical or web
interface, and Phase 6 does not add one.

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

[1] First action
[2] Second action
[0] Back
```

`ui.style` provides `header()`, `section()` and `field()` for these.
Labels in a key/value block are padded to `FIELD_WIDTH` (22) so aligned
blocks line up across modules.

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

### Known deviation

`analyses/mutation.py::_prompt_sequence` offers no `[0]`. It used to mean
"enter a sequence manually", which is a choice rather than an exit; that
moved to `[M]`. A real cancel would need `_prompt_sequence` to gain a
cancel contract and its four call sites to handle it, which is a control
flow change rather than a presentation one. Recorded as follow-up work
rather than done quietly.

## 4. Input

| Helper | Use |
|---|---|
| `ask_choice()` | A menu selection. Renders `Selection: `. |
| `ask(label, default=...)` | A value. A default is shown in brackets and used for empty input. |
| `confirm(question)` | A yes/no question. Anything but an explicit yes means no. |
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

`core/pipeline.py` prints status to the terminal from inside `core/`. It
takes `ui.style` through function-local imports so `core` keeps no
module-level dependency on `ui`. Moving `print_status()` into `ui/` would
be the cleaner arrangement and is recorded as follow-up work; the layering
inversion is confined to two imports rather than spread through the file.

The analysis modules still own their own menus (`run_*_ui`, `*_menu`).
`core/api.py` is the non-interactive layer for every analysis, so a
future interface does not have to go through these menus. Separating the
menus from the analysis modules is the remaining part of audit item P2-2.

## 7. Enforcement

`tests/test_terminal_ui.py` checks the project, not a sample:

1. No separator literal outside `ui/style.py`, with the documented
   `orf.py` table exception.
2. Every `[0]` entry is `Back`, `Cancel` or `Exit`.
3. `Exit` appears only in the main menu and the sequence selector.
4. No `Error:`/`Warning:` prefix variants in user-facing messages.
5. `ui/style.py` imports nothing from `analyses/` or `core/`, so the
   presentation layer stays free of scientific logic.
6. Menu numbers are contiguous, and every entry has a dispatch branch
   and every branch an entry (BLAST menu and main menu).

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
