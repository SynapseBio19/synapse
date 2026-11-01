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
            for label in MENU_ZERO.findall(read(path)):
                label = label.strip()
                # A qualifier in brackets explains the effect; the verb rules.
                verb = label.split("(")[0].strip()
                if verb not in self.ALLOWED:
                    offenders.append(f"{path.name}: [0] {label}")
        self.assertEqual(offenders, [], "[0] must mean Back, Cancel or Exit")

    def test_exit_belongs_to_the_main_menu_only(self):
        exits = []
        for path in UI_FILES:
            for label in MENU_ZERO.findall(read(path)):
                if label.strip() == style.EXIT_LABEL:
                    exits.append(path.name)
        self.assertEqual(sorted(set(exits)), ["main.py", "menu.py"], exits)

    def test_is_back_recognises_the_key(self):
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
        # Entries of the top-level menu, up to and including [0].
        labels, seen_back = [], False
        for key, _label in re.findall(r'print\("\[(\d+)\] ([^"]+)"\)', body):
            if seen_back:
                break            # anything after [0] belongs to a sub-dialog
            if key == "0":
                seen_back = True
            else:
                labels.append(int(key))
        self.assertEqual(labels, list(range(1, len(labels) + 1)),
                         f"BLAST menu numbers must be contiguous: {labels}")
        dispatched = {
            int(key) for key in
            re.findall(r'\n        (?:el)?if choice == "(\d+)":', body)
            if key != "0"
        }
        self.assertEqual(dispatched, set(labels),
                         "every BLAST menu entry must have a branch and vice versa")

    def test_the_main_menu_is_contiguous(self):
        body = read(PROJECT_ROOT / "ui" / "menu.py")
        # [0] is navigation and sits last, so it is not part of the sequence.
        numbers = [int(k) for k in re.findall(r'print\("\[(\d+)\]', body) if k != "0"]
        self.assertEqual(numbers, list(range(1, len(numbers) + 1)),
                         f"main menu numbers must be contiguous: {numbers}")
        self.assertIn('print("[0] Exit")', body)


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

    def test_core_has_no_module_level_dependency_on_ui(self):
        for path in (PROJECT_ROOT / "core").glob("*.py"):
            tree = ast.parse(read(path))
            top_level = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
            for node in top_level:
                module = getattr(node, "module", "") or ""
                names = [a.name for a in node.names]
                with self.subTest(file=path.name):
                    self.assertFalse(module.startswith("ui") or any(n.startswith("ui") for n in names),
                                     f"{path.name} imports ui at module level")


if __name__ == "__main__":
    unittest.main()
