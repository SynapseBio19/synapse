"""The application entry points (Phase 9, 2026-11-19).

The session moved from ``main.py`` to ``synapse_gen/app.py`` so that it can
be installed with the command ``synapse``. ``main.py`` stays as the way to
start Synapse from the source folder; importing it yields the app module
itself, so tests and scripts that patch ``main`` reach the session.
"""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock

PROJECT_ROOT = Path(__file__).resolve().parent.parent


class EntryPointTests(unittest.TestCase):
    def test_main_is_the_app_module(self):
        import main
        from synapse_gen import app
        self.assertIs(main, app)
        self.assertTrue(Path(main.default_sequence_file).is_file())

    def test_cli_ends_cleanly_when_the_input_ends(self):
        from synapse_gen import app
        with mock.patch.object(app, "main", side_effect=EOFError), mock.patch("builtins.print"):
            self.assertEqual(app.cli(), 0)
        with mock.patch.object(app, "main", side_effect=KeyboardInterrupt), mock.patch("builtins.print"):
            self.assertEqual(app.cli(), 130)
        with mock.patch.object(app, "main", return_value=None):
            self.assertEqual(app.cli(), 0)

    def test_python_main_py_quits_without_a_traceback(self):
        for script in ("Q\n", "\n1\n"):          # Q at the file prompt; input ending in the menu
            with self.subTest(script=script):
                completed = subprocess.run([sys.executable, "main.py"], cwd=PROJECT_ROOT, input=script,
                                           capture_output=True, text=True, timeout=120)
                self.assertEqual(completed.returncode, 0, completed.stderr[-500:])
                self.assertNotIn("Traceback", completed.stdout + completed.stderr)

    def test_the_command_is_declared(self):
        text = (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        self.assertIn('synapse = "synapse_gen.app:cli"', text)
        self.assertIn('data = ["examples/*.fasta"]', text)


if __name__ == "__main__":
    unittest.main()
