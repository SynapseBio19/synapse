# Stage F – Independent Closing Review

Date: 2026-11-03. Scope: git history `623a0a8..1516339` (11 commits: the ten
commits of the second Stage F round plus the Enter-key fix `1516339`), all
sixteen interface files, `tests/test_terminal_ui.py`, and the Stage F
documents. No source code was changed for this review. Nothing was sent to
an external service.

Verdict: **Stage F is substantially done, but it is not "fully unified"
yet.** Every interface file really does use the shared renderers. No
module was only added to `SHARED_UI_MODULES`. What remains:

- 59 measured bypasses, none of them documented;
- a core→ui dependency that still exists at runtime;
- stale menu numbers in `FEATURE_STATUS.md`;
- rendered tests that reach less than half of the menus.

A behavior regression (Enter on bracketed prompts) existed from round 1
until `1516339`. Earlier statements that "behavior is preserved" were
therefore wrong for that period.

## 1. Renderer adoption (question 1)

**Confirmed.**

- All 16 files call `ui/style` (about 1,580 calls in total; for example
  `field` 779, `section` 172, `field_line` 142, `header` 102, `ask` 74,
  `menu` 36 plus 2 in `ui/menu.py`).
- No file defines its own rule or header helper. No file prints a message
  prefix or a yes/no notation of its own. Test item 8 enforces both.
- All 11 entries in `SHARED_UI_MODULES` are actually converted. The
  other five files (`main.py`, `ui/report.py`, `core/pipeline.py`,
  `sequence.py`, `cloning.py`) are covered through `INTERFACE_SOURCES`.
- All 25 files changed since `623a0a8` keep their line endings: CRLF
  sources, LF tests and docs. `orf.py` keeps its two intentional bare-LF
  lines.
- `git diff --check` is clean.

**Inconsistent: output that still bypasses the renderers.** The counts
come from an AST scan of literal `print` arguments.

| Kind | Count | Examples |
|---|---|---|
| Validation or status message printed with plain `print` instead of `style.error`/`warning`/`note` | 31 | blast "Invalid database selection.", "Download cancelled.", "No request ID supplied."; cloning "The two enzymes must be different."; crispr "No sequence entered."; literature "Query cannot be empty."; main.py:653 "Please run [2] ORF Analysis first…", main.py:900; restriction "Invalid topology selection." |
| Upper-case heading printed directly instead of `style.section` | 19 | blast "LOCAL DATABASE", "BLAST PROGRAM", "PRIVACY NOTICE", "FASTA RECORDS"; primer "THERMODYNAMIC SETTINGS"; literature "WARNINGS", "SOURCE ERRORS" |
| Dynamic option list built with `print(f"[{i}] …")` plus `style.cancel()` instead of `style.menu` | 9 | blast 3972, 4025, 4350, 4395, 4444, plus the hit lists 4584 and 4611; crispr 2127, 2149, 2404 |
| Field padded by hand | several | rna:2065 "Input source:        Selected ORF"; orf:2994 menu label "Genetic code:        …" |

Enforcement item 8 misses the dynamic lists. Its regex only matches
literal keys such as `[1]`, not f-string keys such as `[{index}]`.

## 2. Core/UI dependency (question 2)

**Not resolved architecturally.** `core/pipeline.py` imports `ui.style`
inside `print_status` (about line 436) and inside the `__main__` block.

- With `ui` blocked through `sys.modules`, importing `core.pipeline` and
  building a pipeline still work, but `print_status()` raises
  `ImportError`.
- The dependency moved from import time to call time. It did not go away.
- `test_core_has_no_module_level_dependency_on_ui` checks only
  module-level imports, so the test confirms the move, not the
  separation.
- `TERMINAL_UI_STANDARD.md` §9 states this openly, which is correct. The
  final report and the change log should not describe it as solved.

Recommendation:

- Move the rendering into `ui/`, for example `ui/report.py:render_status(context)`
  fed by the data `PipelineContext` already exposes.
- Change the only caller (`main.py:888`).
- Add a test that runs `print_status`'s replacement and imports `core`
  with `ui` blocked.

## 3. Variant menu order (question 3)

**Confirmed in code and tests.**

- The new layout is:
  - Input 1–2
  - Processing: 3 Normalize, 4 Filter
  - Analysis: 5 Statistics, 6 Analyze, 7 Compare, 8 Cohort
  - Output 9–11
  - System: 12 Self-Test
- Every dispatch branch was checked against the function it calls:
  - old 8 → new 3;
  - old 6 → new 5;
  - old 3 → new 6;
  - old 5 → new 7;
  - old 7 → new 8.
- No other module, test fixture or golden scenario uses variant menu
  keys.
- `test_variant_menu_follows_the_order_of_the_work` covers the order.

**Inconsistent.**

- `docs/FEATURE_STATUS.md` rows 298–303 still list the old numbers
  (3 Analyze, 5 Compare, 6 Statistics, 7 Cohort, 8 Normalize).
- The same file lists the ORF settings with the numbers from before the
  earlier Stage F reorder: "6.8 Genetic code" and "6.6/6.7 Restore …".
  The menu now has [6] Genetic code, [7] Restore defaults and [8] Restore
  legacy.
- The dispatch branches in `run_variant_analysis_ui` are no longer in
  numeric order (3, 5, 6, 7, 8 are interleaved). This works, but it makes
  checking branch against label harder.

## 4. Behavior preservation (question 4)

The static comparison of all menu entries between `623a0a8` and `HEAD`
shows these differences:

- **From the variant renumbering:** all of them are intended.
- **Entries the static comparison can no longer see:** each now sits in
  a `style.menu` call with a computed list, and each was checked by hand
  and is present:
  - `[F] Load another file` (`ui/menu.py:47`, `zero=None` as before);
  - `[M] Enter sequence manually` (`mutation.py:3691`);
  - topology `[1] Linear`, `[2] Circular` (`restriction.py:2372`);
  - several `[0] Back` (the default `zero`).

Earlier in Stage F, compared with `fc13aad`, there were also these
changes:

- the ORF settings keys moved;
- the Tm menu's `[0]` changed from Back to Cancel.

Both are intentional and documented.

**Intentional behavior changes, all documented:**

- extra pauses after notes or errors in BLAST, variant and CRISPR;
- `confirm` also accepts `yes`/`j`/`ja`;
- the primer Tm dialog asks again on invalid input instead of falling
  back;
- the renumbering of the ORF settings and the variant menu;
- standardized prompt texts ("Selection:" instead of "Select ORF:" /
  "Select topology:").

**Unintended change: the Enter regression, fixed in `1516339`.**

- **Cause:** the automatic rewrite turned `Label [hint]:` prompts into
  `ask(label, default=hint)`, so pressing Enter returned the hint text
  instead of `""`.
- **Prompts affected:** 15 prompts, among them:
  - the literature dates (`none` → "Invalid configuration: date_from …");
  - primer Na⁺, Mg²⁺ and primer concentration;
  - BLAST word size and filter;
  - the cloning strategy count, which looped endlessly.
- **Introduced in:** `5cb4c9a`, `b181ae3`, `fb4070a` and `85c22ac`.
- **Also:** round 1 contained four wrong yes/no conversions
  (`default='y/N'`), fixed in `3853d87`.
- **Status at `HEAD`:** all 15 prompts are restored. The tests in
  `EmptyAnswerTests` fail on the earlier code. A source rule forbids
  `style.ask(..., default=` outside `main.py`.

Remaining risk: other conversions of the same automatic kind may hide
more subtle changes. `ask()` strips the answer, and the prompt texts
changed. The test suite cannot rule this out, because most dialogs are
never driven (see §5).

## 5. Test coverage of rendered screens (question 5)

Measured by instrumenting `style.menu`, `header`, `ask`, `ask_choice` and
`confirm` while the full suite (326 tests) ran:

| Call | Call sites | Executed by any test |
|---|---:|---:|
| `style.menu` | 38 | 17 |
| `style.header` | 102 | 20 |
| `style.ask_choice` | 33 | 16 |
| `style.ask` | 72 | 16 |
| `style.confirm` | 7 | 1 |

These 21 menus are never rendered by any test:

| File | Lines |
|---|---|
| blast | 3997, 4094, 4289, 4380 |
| cloning | 4407, 4565, 4669, 4710 |
| crispr | 2317, 2367, 2434, 2459, 2481 |
| mutation | 3688, 3820 |
| orf | 2931, 2987 (the settings menu) |
| primer | 3837 |
| main.py | 225, 262 |
| ui/menu.py | 47 (the sequence selector) |

What is covered:

- the first screen of the 11 module menus and the main menu (frame,
  title, context fields, footer, unknown key);
- two report screens;
- the confirm notation;
- three empty-answer dialogs.

What is not covered:

- most submenus;
- result views and tables;
- error screens such as invalid FASTA or an empty query;
- cancellation paths;
- the nine dynamic selection lists.

Several of these are tested only through source rules. A source rule
checks how the code is written, not what the screen shows. The coverage
is therefore **not sufficient** to support the claim "behavior preserved
in all dialogs".

## 6. Documented exceptions (question 6)

`TERMINAL_UI_STANDARD.md` §9 lists seven exceptions. Each was checked.

- **Raw line input:** necessary. `crispr:2169` and `variant:2472` / `2616`
  read pasted lines unchanged. Not complete: `crispr:2221`
  (`_read_local_references`) also uses `input()`, but it strips the line,
  so the stated reason ("must not be stripped") does not apply to it.
  It could use `style.ask`, or §9 should name it.
- **Two CRISPR numeric helpers** (`crispr:2098`, `2113`): accurate.
  Whether they are necessary is debatable: changing the callers' prompt
  strings is a mechanical change.
- **Closing rules:** the count is accurate. There are 18 direct
  `print(major_rule())` calls (main 1, report 1, pipeline 1, orf 2,
  primer 7, rna 3, protein 2, restriction 1).
- **Result lines:** necessary and appropriate.
- **Mutation `zero=None`:** accurate. The same applies to the sequence
  selector, which §9 does not mention. That selector is listed in §3
  instead.
- **`core/pipeline.py`:** accurate (see §2).
- **Menus inside analysis modules (P2-2):** accurate.

**Not documented anywhere:**

- the 31 plain messages;
- the 19 direct headings;
- the 9 dynamic option lists;
- the hand-padded fields.

None of them has a reason that needs an exception. They are leftovers of
the conversion, so the list of exceptions is **not complete**.

## 7. Reproducibility (question 7)

The environment:

| | |
|---|---|
| OS | Windows 11 |
| Python | 3.14.6, `.venv` |
| pyright | 1.1.414, run as `python -m pyright --pythonversion 3.10`. The `pyright.exe` script is blocked by application control. |

| Check | Result |
|---|---|
| `python -m unittest discover -s tests -t .`, run 1 | 326 tests, OK, 34.9 s |
| same, run 2 | 326 tests, OK, 35.2 s |
| Golden master (`tests.test_golden_master`) | 12 tests, OK |
| Golden reference files changed since `623a0a8` | none (`git log 623a0a8..HEAD -- tests/golden tests/golden_scenarios.py` is empty) |
| Built-in self-tests | `blast.run_self_test`, `crispr.run_self_test`, `literature.self_test`: all `True` |
| pyright | 112 errors; baseline 134; 0 new diagnostics (compared by file, rule and message, ignoring line numbers) |
| `git diff --check 623a0a8 HEAD` (cr-at-eol) | clean |
| Working tree after all runs | clean |

All the figures in the reports could be reproduced. The suite runs no
network test. The figures do not show that the external services work.

## 8. Recommendations, in order

1. **Correct the claims.**
   - Do not call Stage F "fully unified" until items 2–4 below are done.
   - State in `PHASE6_FINAL_REPORT.md` and `CHANGES.md` that the
     behavior of round 1 and round 2 was not preserved until `1516339`.
2. **Convert the 31 messages, 19 headings and 9 dynamic lists** to
   `style.error`/`warning`/`note`, `style.section` and `style.menu`.
   - Extend enforcement item 8 to catch `print(f"[{…}] …")` and plain
     messages that start with "Invalid", "Please", "No ", "Cannot",
     "Unknown" and similar words.
3. **Update the stale numbers** in `FEATURE_STATUS.md`: the variant rows
   298–303 and the ORF settings rows 115–116.
   - Sort the variant dispatch branches numerically.
4. **Remove the core→ui dependency at runtime** as described in §2, with
   a test that imports `core` while `ui` is blocked.
5. **Extend the rendered tests** to every `style.menu` call site, using
   the instrumentation from §5 as a coverage check (38 of 38).
   - Add at least one test per module that drives an error screen, a
     cancellation and an invalid entry.
   - Drive every `confirm` site and the dialogs that read numbers with
     an empty answer.
6. **Settle §9:** add `crispr:2221` and the sequence selector, or remove
   the exceptions by converting the code.
