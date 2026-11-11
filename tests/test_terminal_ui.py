"""Stage F: the terminal interface standard is enforced, not just documented.

Every check scans the whole project rather than a sample, so a module that
drifts from docs/TERMINAL_UI_STANDARD.md fails the suite.
"""

from __future__ import annotations

import ast
import re
import unittest
from pathlib import Path

from tests import _support  # noqa: F401
from tests._support import PROJECT_ROOT

from ui import style

# Every file that draws a terminal interface.
UI_FILES = sorted(
    [p for p in (PROJECT_ROOT / "analyses").glob("*.py") if p.name != "__init__.py"]
    + [p for p in (PROJECT_ROOT / "ui").glob("*.py") if p.name != "__init__.py"]
    + [PROJECT_ROOT / "core" / "pipeline.py", PROJECT_ROOT / "main.py"]
)

SEPARATOR = re.compile(r'"(=|-)"\s*\*\s*')
MENU_ZERO = re.compile(r'print\(\s*f?"\s*\[0\]\s*([^"]*)"')
BAD_PREFIX = re.compile(r'print\(\s*f?"(Error|Warning|Note)[: ]')


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def zero_labels(source: str) -> list:
    """Every label given to [0] in ``source``: printed or through ui.style."""
    labels = [label.strip() for label in MENU_ZERO.findall(source)]
    constants = {"BACK_LABEL": style.BACK_LABEL, "CANCEL_LABEL": style.CANCEL_LABEL,
                 "EXIT_LABEL": style.EXIT_LABEL}

    def text(node, default):
        if node is None:
            return default
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Attribute) and node.attr in constants:
            return constants[node.attr]
        return "<expression>"

    for node in ast.walk(ast.parse(source)):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "style"):
            continue
        name = node.func.attr
        if name == "menu":
            zero = next((k.value for k in node.keywords if k.arg == "zero"), None)
            has_zero = any(k.arg == "zero" for k in node.keywords)
            # A conditional [0] (present or absent) is checked on both branches.
            branches = [zero.body, zero.orelse] if isinstance(zero, ast.IfExp) else [zero]
            for branch in branches:
                label = text(branch, style.BACK_LABEL) if has_zero else style.BACK_LABEL
                if label is not None:
                    labels.append(label)
        elif name == "back":
            labels.append(text(node.args[0] if node.args else None, style.BACK_LABEL))
        elif name == "cancel":
            labels.append(style.CANCEL_LABEL)
        elif name == "exit_option":
            labels.append(style.EXIT_LABEL)
        elif name == "option" and node.args and isinstance(node.args[0], ast.Constant) \
                and str(node.args[0].value) == "0":
            labels.append(text(node.args[1], ""))
    return labels


def first_menu_keys(source: str) -> list:
    """The entry keys of the first style.menu([...]) call in ``source``."""
    for node in ast.walk(ast.parse(source)):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "menu" and node.args):
            groups = ast.literal_eval(node.args[0])
            return [str(key) for _title, entries in groups for key, _label in entries]
    return []


class SeparatorTests(unittest.TestCase):
    def test_no_module_draws_its_own_separator(self):
        offenders = []
        for path in UI_FILES:
            if path.name == "style.py":
                continue
            for number, line in enumerate(read(path).splitlines(), 1):
                if SEPARATOR.search(line):
                    # Documented exception: a table rule spans its own table.
                    if path.name == "orf.py" and "24 + 14 * len(candidates)" in line:
                        continue
                    offenders.append(f"{path.name}:{number}: {line.strip()}")
        self.assertEqual(offenders, [], "separator literals must come from ui.style")

    def test_one_width_for_the_whole_application(self):
        self.assertEqual(style.RULE_WIDTH, 70)
        self.assertEqual(len(style.major_rule()), style.RULE_WIDTH)
        self.assertEqual(len(style.minor_rule()), style.RULE_WIDTH)
        self.assertNotEqual(style.major_rule(), style.minor_rule())

    def test_every_ui_file_uses_the_shared_rules(self):
        users = [p.name for p in UI_FILES
                 if "major_rule()" in read(p) or "minor_rule()" in read(p)]
        # Every module that frames a screen must go through ui.style.
        self.assertGreaterEqual(len(users), 15, users)


class NavigationTests(unittest.TestCase):
    ALLOWED = {style.BACK_LABEL, style.CANCEL_LABEL, style.EXIT_LABEL}

    def test_zero_always_leaves_the_screen(self):
        offenders = []
        for path in UI_FILES:
            for label in zero_labels(read(path)):
                # A qualifier in brackets explains the effect; the verb rules.
                verb = label.split("(")[0].strip()
                if verb not in self.ALLOWED:
                    offenders.append(f"{path.name}: [0] {label}")
        self.assertEqual(offenders, [], "[0] must mean Back, Cancel or Exit")

    def test_exit_belongs_to_the_main_menu_only(self):
        exits = []
        for path in UI_FILES:
            for label in zero_labels(read(path)):
                if label == style.EXIT_LABEL:
                    exits.append(path.name)
        self.assertEqual(sorted(set(exits)), ["main.py", "menu.py"], exits)

    def test_is_back_recognizes_the_key(self):
        self.assertTrue(style.is_back("0"))
        self.assertTrue(style.is_back("  0 "))
        self.assertFalse(style.is_back("1"))
        self.assertFalse(style.is_back(""))


class MenuNumberingTests(unittest.TestCase):
    """A menu's numbers must run in order and every one must be handled."""

    def _menu_source(self, module: str, function: str) -> str:
        source = read(PROJECT_ROOT / module)
        tree = ast.parse(source)
        node = next(n for n in tree.body
                    if isinstance(n, ast.FunctionDef) and n.name == function)
        return "\n".join(source.split("\n")[node.lineno - 1:node.end_lineno])

    def test_blast_menu_is_contiguous_and_fully_dispatched(self):
        body = self._menu_source("analyses/blast.py", "blast_menu")
        # Entries of the top-level menu: the first style.menu() call.
        labels = [int(key) for key in first_menu_keys(body)]
        self.assertEqual(labels, list(range(1, len(labels) + 1)),
                         f"BLAST menu numbers must be contiguous: {labels}")
        dispatched = {
            int(key) for key in
            re.findall(r'\n        (?:el)?if choice == "(\d+)":', body)
            if key != "0"
        }
        self.assertEqual(dispatched, set(labels),
                         "every BLAST menu entry must have a branch and vice versa")

    def test_variant_menu_follows_the_order_of_the_work(self):
        """Master Prompt 5.3: input, normalization, filtering, statistics, interpretation."""
        body = self._menu_source("analyses/variant_analysis.py", "run_variant_analysis_ui")
        menu = next(n for n in ast.walk(ast.parse(body)) if isinstance(n, ast.Call)
                    and isinstance(n.func, ast.Attribute) and n.func.attr == "menu")
        labels = {label: int(key) for _title, entries in ast.literal_eval(menu.args[0])
                  for key, label in entries}
        keys = sorted(labels.values())
        self.assertEqual(keys, list(range(1, len(keys) + 1)))
        order = ["Load VCF", "Normalize Variants", "Filter Variants", "Variant Statistics",
                 "Analyze Variants", "Full Variant Report"]
        self.assertEqual([labels[name] for name in order], sorted(labels[name] for name in order))
        dispatched = {int(k) for k in re.findall(r'choice == "(\d+)"', body) if k != "0"}
        self.assertEqual(dispatched, set(keys))

    def test_the_main_menu_is_contiguous(self):
        body = read(PROJECT_ROOT / "ui" / "menu.py")
        # Letter keys (P, R, S) are actions; [0] is navigation and sits last.
        numbers = [int(k) for k in first_menu_keys(body) if k.isdigit()]
        self.assertEqual(numbers, list(range(1, len(numbers) + 1)),
                         f"main menu numbers must be contiguous: {numbers}")
        self.assertIn("zero=style.EXIT_LABEL", body)


class MessageTests(unittest.TestCase):
    def test_message_prefixes_are_uniform(self):
        offenders = []
        for path in UI_FILES:
            for number, line in enumerate(read(path).splitlines(), 1):
                if BAD_PREFIX.search(line):
                    offenders.append(f"{path.name}:{number}: {line.strip()}")
        self.assertEqual(offenders, [], "use ERROR:/WARNING:/NOTE: (see ui.style)")

    def test_helpers_render_the_documented_prefixes(self):
        import contextlib
        import io
        for helper, prefix in ((style.error, "ERROR:"), (style.warning, "WARNING:"),
                               (style.note, "NOTE:"), (style.success, "OK:")):
            with self.subTest(prefix=prefix):
                buffer = io.StringIO()
                with contextlib.redirect_stdout(buffer):
                    helper("something happened")
                self.assertEqual(buffer.getvalue().strip(), f"{prefix} something happened")


class LayeringTests(unittest.TestCase):
    def test_style_module_carries_no_scientific_logic(self):
        tree = ast.parse(read(PROJECT_ROOT / "ui" / "style.py"))
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        self.assertEqual(imported & {"analyses", "core", "data", "fileio"}, set(),
                         f"ui.style must stay presentation-only, imports: {sorted(imported)}")

    def test_core_imports_ui_nowhere(self):
        """Not at module level and not inside a function: a deferred import is still a dependency."""
        hits = []
        for path in sorted((PROJECT_ROOT / "core").glob("*.py")):
            for node in ast.walk(ast.parse(read(path))):
                if isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                    modules = [node.module]
                elif isinstance(node, ast.Import):
                    modules = [alias.name for alias in node.names]
                else:
                    continue
                hits.extend(f"{path.name}:{node.lineno}: {m}" for m in modules
                            if m == "ui" or m.startswith("ui."))
        self.assertEqual(hits, [])

    def test_core_runs_with_the_ui_package_unavailable(self):
        """Import every core module and run core's own logic while `ui` cannot be imported."""
        import subprocess
        import sys
        import textwrap
        probe = textwrap.dedent('''
            import importlib.abc, json, pkgutil, sys, tempfile, os

            class BlockUi(importlib.abc.MetaPathFinder):
                def find_spec(self, name, path, target=None):
                    if name == "ui" or name.startswith("ui."):
                        raise ImportError("ui is unavailable: " + name)

            sys.meta_path.insert(0, BlockUi())
            import core
            for info in pkgutil.iter_modules(core.__path__):
                __import__("core." + info.name)

            from core.pipeline import SynapsePipeline
            from core.report import build_analysis_report, export_analysis_report
            from core.alphabets import classify_sequence
            from core.seqtools import reverse_complement, get_genetic_code

            pipeline = SynapsePipeline("ATGAAACCCGGGTTTTAA" * 5, "demo")
            report = pipeline.status_report()
            assert report["length"] == 90 and len(report["analyses"]) == 12, report
            seen = []
            pipeline.print_status(seen.append)          # any renderer, here a list
            assert seen == [report]
            assert classify_sequence("ACGU").sequence_type == "RNA"
            assert reverse_complement("AACG") == "CGTT"
            assert get_genetic_code(11).translate("ATGAAATAA") == "MK"
            build_analysis_report(pipeline)
            export_analysis_report(pipeline, os.path.join(tempfile.mkdtemp(), "r.json"))
            assert not any(n == "ui" or n.startswith("ui.") for n in sys.modules)
            print("CORE WITHOUT UI OK")
        ''')
        completed = subprocess.run([sys.executable, "-c", probe], cwd=PROJECT_ROOT,
                                   capture_output=True, text=True, timeout=120)
        self.assertIn("CORE WITHOUT UI OK", completed.stdout, completed.stderr[-2000:])

    def test_core_smoke_test_runs_without_ui(self):
        """`python -m core.pipeline` prints the status as data, with ui unavailable."""
        import json
        import subprocess
        import sys
        probe = ("import importlib.abc, runpy, sys\n"
                 "class B(importlib.abc.MetaPathFinder):\n"
                 "    def find_spec(self, name, path, target=None):\n"
                 "        if name == 'ui' or name.startswith('ui.'): raise ImportError(name)\n"
                 "sys.meta_path.insert(0, B())\n"
                 "runpy.run_module('core.pipeline', run_name='__main__')\n")
        completed = subprocess.run([sys.executable, "-c", probe], cwd=PROJECT_ROOT,
                                   capture_output=True, text=True, timeout=120)
        self.assertEqual(completed.returncode, 0, completed.stderr[-2000:])
        self.assertEqual(json.loads(completed.stdout)["header"], "TEST_DNA")

    def test_terminal_status_screen_is_drawn_by_ui(self):
        """The renderer core receives is ui.report's, and it draws the shared screen."""
        import contextlib
        import io
        from core.pipeline import SynapsePipeline
        from ui.report import print_pipeline_status
        pipeline = SynapsePipeline("ACGT" * 30, "demo")
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            pipeline.print_status(print_pipeline_status)
        lines = buffer.getvalue().splitlines()
        self.assertEqual(lines[1:4], [style.major_rule(), "PIPELINE STATUS", style.major_rule()])
        self.assertIn(style.field_line("Length", "120 bp"), lines)
        self.assertIn(style.columns(("Sequence Overview", "NOT RUN"), (28,)), lines)
        self.assertEqual(lines[-1], style.major_rule())


# ----------------------------------------------------------------------
# The rendered interface, not only source patterns (Stage F completion)
# ----------------------------------------------------------------------
# Menu screens rendered by the checks below: every module menu and the
# main menu. Each entry is converted to ui.style, not only listed here.
SHARED_UI_MODULES = ("orf", "crispr", "blast", "literature", "primer", "variant_analysis",
                     "rna", "protein", "restriction", "mutation", "menu")

# Every file that draws terminal screens; all of them are routed through
# ui.style, so the source rule applies to all.
INTERFACE_SOURCES = tuple(
    [PROJECT_ROOT / "main.py", PROJECT_ROOT / "core" / "pipeline.py"]
    + [PROJECT_ROOT / "ui" / name for name in ("menu.py", "report.py")]
    + [PROJECT_ROOT / "analyses" / f"{name}.py" for name in (
        "orf", "crispr", "blast", "literature", "primer", "variant_analysis", "sequence",
        "rna", "protein", "restriction", "cloning", "mutation")]
)

# Answers that reach the menu: RNA first shows its pipeline input and
# waits for Enter.
LEAD_IN = {"rna": ("",)}


def render_first_screen(module, answers=("0",)):
    """The text a module prints until its main menu has been left."""
    import builtins
    import contextlib
    import io
    import logging
    from unittest import mock
    from tests._support import deterministic_dna

    dna = deterministic_dna(900)
    pending = list(answers)

    def scripted(prompt=""):
        print(prompt, end="")
        if not pending:
            raise EOFError("script exhausted")
        answer = pending.pop(0)
        print(answer)
        return answer

    calls = {
        "orf": lambda: __import__("analyses.orf", fromlist=["x"]).run_orf_analysis_ui(dna),
        "crispr": lambda: __import__("analyses.crispr", fromlist=["x"]).crispr_menu(dna, "demo"),
        "blast": lambda: __import__("analyses.blast", fromlist=["x"]).blast_menu(dna, "demo"),
        "literature": lambda: __import__("analyses.literature", fromlist=["x"])._interactive_menu(),
        "primer": lambda: __import__("analyses.primer", fromlist=["x"]).run_primer_analysis_ui(dna, "demo"),
        "variant_analysis": lambda: __import__("analyses.variant_analysis", fromlist=["x"]).run_variant_analysis_ui(),
        "rna": lambda: __import__("analyses.rna", fromlist=["x"]).run_rna_analysis_ui(dna, header="demo"),
        "protein": lambda: __import__("analyses.protein", fromlist=["x"]).run_protein_analysis_from_record(
            "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAAEERNWLGAMLGRSAL",
            name="demo"),
        "restriction": lambda: __import__("analyses.restriction", fromlist=["x"]).run_restriction_cloning_ui(
            dna, "demo", [{"header": "demo", "sequence": dna, "topology": "linear"}]),
        "mutation": lambda: __import__("analyses.mutation", fromlist=["x"]).run_mutation_variant_ui(
            dna, "demo", [{"header": "demo", "sequence": dna}]),
        "menu": lambda: __import__("ui.menu", fromlist=["x"]).show_menu(),
    }
    buffer = io.StringIO()
    logging.disable(logging.CRITICAL)
    try:
        with mock.patch.object(builtins, "input", scripted), contextlib.redirect_stdout(buffer):
            try:
                calls[module]()
            except EOFError:
                # The script ran out before the screen was left; the
                # caller's assertions then see where it stopped.
                print("<<script exhausted>>")
    finally:
        logging.disable(logging.NOTSET)
    return buffer.getvalue()


class RenderedScreenTests(unittest.TestCase):
    # A module's menu carries the name of its main-menu entry.
    TITLES = {"orf": "ORF ANALYSIS", "crispr": "CRISPR GUIDE DESIGN", "blast": "BLAST SEARCH",
              "literature": "LITERATURE SEARCH", "primer": "PRIMER DESIGN",
              "variant_analysis": "VARIANT ANALYSIS", "rna": "RNA ANALYSIS",
              "protein": "PROTEIN ANALYSIS", "restriction": "RESTRICTION & CLONING ANALYSIS",
              "mutation": "MUTATION & VARIANT ANALYSIS", "menu": "MAIN MENU"}
    ZERO = {"menu": "[0] Exit"}

    def _first_screen(self, module, answer="0"):
        lines = render_first_screen(module, LEAD_IN.get(module, ()) + (answer,)).splitlines()
        end = next(i for i, line in enumerate(lines) if line.startswith(style.SELECTION_PROMPT + ":")
                   and line.endswith(answer))
        # The menu screen starts at the last title frame before its prompt.
        start = max(i for i in range(end - 2) if lines[i] == style.major_rule()
                    and lines[i + 2] == style.major_rule())
        return lines[start:end + 1]

    def test_each_screen_has_the_shared_frame_title_and_footer(self):
        for module in SHARED_UI_MODULES:
            with self.subTest(module=module):
                screen = self._first_screen(module)
                self.assertEqual(screen[:3], [style.major_rule(), self.TITLES[module], style.major_rule()])
                self.assertEqual(screen[-1], "Selection: 0")
                self.assertEqual(screen[-2], style.minor_rule())
                self.assertEqual(screen[-3], self.ZERO.get(module, "[0] Back"))
                self.assertEqual(screen[-4], "")

    def test_context_fields_line_up(self):
        field = re.compile(r"^[A-Z][^:\[]*:\s")
        for module in SHARED_UI_MODULES:
            with self.subTest(module=module):
                screen = self._first_screen(module)
                first_option = next(i for i, line in enumerate(screen) if line.startswith("["))
                context = [line for line in screen[3:first_option] if field.match(line)]
                for line in context:
                    label, _, _ = line.partition(":")
                    if len(label) + 1 <= style.FIELD_WIDTH:
                        self.assertRegex(line, rf"^.{{{style.FIELD_WIDTH}}} \S", line)

    def test_an_unknown_key_gets_the_shared_reply(self):
        # The main menu returns its answer to main.py, which replies.
        for module in [m for m in SHARED_UI_MODULES if m != "menu"]:
            with self.subTest(module=module):
                output = render_first_screen(module, LEAD_IN.get(module, ()) + ("x", "0"))
                after = output.split("Selection: x\n", 1)[1]
                self.assertEqual(after.lstrip("\n").splitlines()[0], style.INVALID_SELECTION)

    def test_converted_modules_do_not_bypass_ui_style(self):
        """No message prefix, menu entry, selection prompt or rule helper of their own."""
        rules = {
            "message prefix printed directly": re.compile(r'print\(f?"(\\n)*(ERROR|WARNING|NOTE|OK): '),
            "menu entry printed directly": re.compile(r'print\(f?"(\\n)*\[[0-9A-Z]+\] '),
            "selection prompt of its own": re.compile(r'input\(f?"(\\n)*(Selection|Select option)'),
            "own rule or header helper": re.compile(r"def (_line|_print_\w*header|_\w+_ui_header)\("),
            "SYNAPSE title prefix": re.compile(r'header\(f?["\']SYNAPSE'),
        }
        rules["yes/no question of its own"] = re.compile(r'\[[yY]/[nN]\]|\([yY]/[nN]\)')
        hits = []
        for path in INTERFACE_SOURCES:
            for number, line in enumerate(read(path).splitlines(), 1):
                hits.extend(f"{path.name}:{number}: {name}" for name, rule in rules.items()
                            if rule.search(line))
        self.assertEqual(hits, [])

    def test_report_screens_align_their_fields(self):
        """The sequence overview and the pipeline status, as rendered."""
        import contextlib
        import io
        import os
        import subprocess
        import sys
        import tempfile
        from core.pipeline import SynapsePipeline
        from ui.report import print_pipeline_status
        with tempfile.TemporaryDirectory() as logs:
            env = dict(os.environ, SYNAPSE_LOG_DIR=logs, PYTHONIOENCODING="utf-8")
            overview = subprocess.run([sys.executable, "main.py"], cwd=PROJECT_ROOT, input="\n1\n0\n",
                                      capture_output=True, text=True, encoding="utf-8",
                                      timeout=120, env=env).stdout
        buffer = io.StringIO()
        pipeline = SynapsePipeline("ACGT" * 30, "demo")
        with contextlib.redirect_stdout(buffer):
            pipeline.print_status(print_pipeline_status)
        field = re.compile(r"^([A-Za-z][^:\[]{0,20}):( +)\S")
        for name, text in (("sequence overview", overview), ("pipeline status", buffer.getvalue())):
            with self.subTest(screen=name):
                lines = text.splitlines()
                self.assertIn(style.major_rule(), lines)
                fields = [line for line in lines if field.match(line)
                          and not line.startswith(("ERROR:", "WARNING:", "NOTE:", "OK:"))]
                self.assertGreater(len(fields), 3)
                for line in fields:
                    self.assertRegex(line, rf"^.{{{style.FIELD_WIDTH}}} \S", line)
                headings = [lines[i - 1] for i, line in enumerate(lines)
                            if line == style.minor_rule() and i and lines[i - 1].strip()
                            and not lines[i - 1].startswith(("[", "Selection"))]
                for heading in headings:
                    self.assertEqual(heading, heading.upper(), heading)

    def test_yes_no_questions_use_one_notation(self):
        import builtins
        import contextlib
        import io
        from unittest import mock
        for answer, default_yes, expected in (("", False, False), ("y", False, True), ("yes", False, True),
                                              ("", True, True), ("n", True, False)):
            buffer = io.StringIO()
            with mock.patch.object(builtins, "input", lambda prompt: (print(prompt, end=""), answer)[1]), \
                    contextlib.redirect_stdout(buffer):
                result = style.confirm("Proceed?", default_yes=default_yes)
            self.assertEqual(result, expected, (answer, default_yes))
            self.assertEqual(buffer.getvalue(), "Proceed? [Y/n]: " if default_yes else "Proceed? [y/N]: ")

    def test_menu_renderer_layout(self):
        import contextlib
        import io
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            style.menu([("Search", [(1, "Run")]), (None, [(2, "Other")])], zero="Cancel")
        self.assertEqual(buffer.getvalue().splitlines(),
                         ["", "SEARCH", "[1] Run", "", "[2] Other", "", "[0] Cancel", style.minor_rule()])
        self.assertEqual(style.columns(["ORF-001", "+3", "87 bp"], [8, 6]), "ORF-001  +3     87 bp")
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            style.menu([(None, [("F", "Load another file")])], zero=None)
        self.assertEqual(buffer.getvalue().splitlines(), ["", "[F] Load another file", "", style.minor_rule()])



# ----------------------------------------------------------------------
# Pressing Enter keeps the documented default (regression, 2026-11-03)
# ----------------------------------------------------------------------
class _Stop(Exception):
    """Raised by a stub to end a dialog once it has handed over its input."""


def _scripted(answers):
    pending = list(answers)

    def answer(prompt=""):
        if not pending:
            raise EOFError("script exhausted")
        return pending.pop(0)
    return answer


class EmptyAnswerTests(unittest.TestCase):
    """The Stage F conversion once turned bracket hints into values:
    Enter at "From date ... [none]" returned "none" and every literature
    search failed with "Invalid configuration"; "[50.0 mM]" broke the primer
    conditions and "[default 5]" the cloning planner. Enter must give the
    caller an empty answer, which each caller maps to its own default."""

    def _run(self, answers, call):
        import builtins
        import contextlib
        import io
        from unittest import mock
        with mock.patch.object(builtins, "input", _scripted(answers)), \
                contextlib.redirect_stdout(io.StringIO()):
            call()

    def test_literature_search_accepts_enter_for_every_option(self):
        from unittest import mock
        from analyses import literature
        seen = []

        def stub(query):
            seen.append(query)
            raise _Stop

        with mock.patch.object(literature, "run_literature_search", stub):
            with self.assertRaises(_Stop):
                self._run(["1", "BRCA1 DNA repair", "", "", "", "", "", ""],
                          literature._interactive_menu)
        query = seen[0]
        self.assertEqual((query.date_from, query.date_to, query.max_results),
                         (None, None, 30))
        self.assertFalse(query.open_access_only)
        self.assertEqual(query.sort_by, "relevance")

    def test_primer_conditions_accept_enter_for_the_defaults(self):
        from analyses import primer
        result = []
        self._run(["", "", "", ""], lambda: result.append(primer._primer_ui_read_conditions()))
        conditions = result[0]
        self.assertIsNotNone(conditions)
        self.assertEqual((conditions.sodium_mM, conditions.magnesium_mM, conditions.dntp_mM,
                          conditions.primer_nM), (50.0, 1.5, 0.2, 250.0))

    def test_cloning_planner_accepts_enter_for_the_strategy_count(self):
        from unittest import mock
        from analyses import cloning
        seen = []

        def stub(*args, **kwargs):
            seen.append(kwargs.get("top_n"))
            raise _Stop

        with mock.patch.object(cloning, "automatic_cloning_planner", stub):
            with self.assertRaises(_Stop):
                self._run(["2", ""], lambda: cloning._cloning_ui_automatic_planner("ACGT" * 50, "ACGT" * 20))
        self.assertEqual(seen, [5])

    def test_primer_tm_model_accepts_enter_for_nearest_neighbor(self):
        """"Selection [1]" took Enter as model 1 before Stage F; the rewrite had lost that key."""
        from analyses import primer
        result = []
        self._run([""], lambda: result.append(primer._primer_ui_select_tm_model()))
        self.assertEqual(result, ["nearest_neighbor"])

    def test_modules_keep_bracket_hints_in_the_label(self):
        """A bracket in a prompt is a hint; ask(..., default=...) would return it as a value."""
        hits = [f"{path.name}:{number}" for path in INTERFACE_SOURCES if path.name != "main.py"
                for number, line in enumerate(read(path).splitlines(), 1)
                if re.search(r"style\.ask\([^)]*default=", line)]
        self.assertEqual(hits, [])


def _leading_text(node):
    """The literal text a print argument starts with, and the node that follows it.

    Returns (text, next_node): for "abc" -> ("abc", None); for f"[{i}] x" ->
    ("[", <i>); for "Unknown: " + x -> ("Unknown: ", None). None if the
    argument does not start with literal text.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value, None
    if isinstance(node, ast.JoinedStr):
        text = ""
        for part in node.values:
            if isinstance(part, ast.FormattedValue):
                return text, part.value
            if isinstance(part, ast.Constant):
                text += str(part.value)
        return text, None
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return _leading_text(node.left)
    return None, None


def direct_output(source: str):
    """(line, kind, text) for each print() that bypasses a shared renderer.

    - "computed menu key": an entry "[{i}] ..." printed inside a loop over
      enumerate(), whose counter i is the key (a hand-built option list);
    - "uncategorized message": a validation or status sentence without
      error()/warning()/note()/success()/invalid_selection();
    - "direct heading": an upper-case heading without section();
    - "rule inside a printed text": a heading drawn with its own dashes.
    """
    tree = ast.parse(source)
    counters = {}
    for loop in ast.walk(tree):
        if isinstance(loop, (ast.For, ast.comprehension)) and isinstance(loop.iter, ast.Call) \
                and getattr(loop.iter.func, "id", None) == "enumerate" \
                and isinstance(loop.target, ast.Tuple) and isinstance(loop.target.elts[0], ast.Name):
            for inner in ast.walk(loop) if isinstance(loop, ast.For) else ():
                counters[id(inner)] = loop.target.elts[0].id
    found = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and getattr(node.func, "id", None) == "print" and node.args):
            continue
        text, following = _leading_text(node.args[0])
        if text is None:
            continue
        stripped = text.lstrip("\r\n ")
        if re.search(r"\n[-=]{3,}", text):
            found.append((node.lineno, "rule inside a printed text", text))
        elif stripped == "[" and isinstance(following, ast.Name) and counters.get(id(node)) == following.id:
            found.append((node.lineno, "computed menu key", text))
        elif DIRECT_MESSAGE.match(stripped) and not any(
                stripped.startswith(allowed) for allowed in RESULT_LINES):
            found.append((node.lineno, "uncategorized message", text))
        elif following is None and DIRECT_HEADING.fullmatch(stripped) \
                or following is not None and stripped.endswith("(") and DIRECT_HEADING.fullmatch(stripped[:-1].strip()):
            found.append((node.lineno, "direct heading", text))
    return sorted(found)


# Sentences that are interface messages, whatever their wording in detail.
DIRECT_MESSAGE = re.compile(
    r"(Invalid |Please |Cannot |Could not |Unknown |Expected |Query cannot|Select at least"
    r"|No \w+(?: \w+)? (?:supplied|entered|selected|loaded)\b"
    r"|IMPORTANT:|Returning to |Resumed |Removed |Update completed"
    r"|[A-Z][\w -]* (?:cancelled|completed|completed successfully)\.$"
    r"|No local specificity analysis|No compatible local)"
)
# Lines that look like messages but are part of a result.
RESULT_LINES = (
    "Invalid character '",  # main.py: the positions listed under the validation ERROR
)
DIRECT_HEADING = re.compile(r"[A-Z][A-Z0-9]+(?: [A-Z0-9/&+-]+)+|[A-Z]{4,}")


class DirectOutputTests(unittest.TestCase):
    """Output the earlier source rule missed: computed keys, plain messages, plain headings."""

    def test_rules_recognize_each_kind(self):
        sample = (
            "for i, name in enumerate(names, 1):\n"
            "    print(f'[{i}] {name}')\n"
            "for result in results:\n"
            "    print(f\"[{'PASS' if ok else 'FAIL'}] {result}\")\n"
            "print('Invalid engine.')\n"
            "print('\\nDownload cancelled.')\n"
            "print('No FASTA path supplied.')\n"
            "print(\"Invalid character '\" + c)\n"
            "print('\\nPRIVACY NOTICE')\n"
            "print(f'RANKED GUIDES ({mode})')\n"
            "print(f'HSP {n}: identity')\n"
            "print('No restriction sites detected.')\n"
            "print('\\nABSTRACT\\n--------')\n"
        )
        self.assertEqual([(line, kind) for line, kind, _ in direct_output(sample)], [
            (2, "computed menu key"), (5, "uncategorized message"), (6, "uncategorized message"),
            (7, "uncategorized message"), (9, "direct heading"), (10, "direct heading"),
            (13, "rule inside a printed text"),
        ])

    def test_no_message_is_drawn_as_a_field(self):
        """style.field('Error', ...) looks like data; an error is an ERROR line."""
        label = re.compile(r"(Error|Warning|Note|.+ (failed|error))")
        hits = []
        for path in INTERFACE_SOURCES:
            for node in ast.walk(ast.parse(read(path))):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                        and node.func.attr == "field" and node.args \
                        and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str) \
                        and label.fullmatch(node.args[0].value):
                    hits.append(f"{path.relative_to(PROJECT_ROOT).as_posix()}:{node.lineno}: {node.args[0].value}")
        self.assertEqual(hits, [])

    def test_no_interface_file_prints_around_the_renderers(self):
        hits = [f"{path.relative_to(PROJECT_ROOT).as_posix()}:{line}: {kind}: {text.strip()!r}"
                for path in INTERFACE_SOURCES for line, kind, text in direct_output(read(path))]
        self.assertEqual(hits, [])


if __name__ == "__main__":
    unittest.main()
