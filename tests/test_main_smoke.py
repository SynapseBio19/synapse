"""End-to-end smoke test: drive main.py through its real menus.

The script runs in a subprocess with scripted keyboard input:
load the example file -> sequence 1 -> ORF analysis, select ORF 1 ->
RNA analysis -> protein analysis -> pipeline status -> exit.
It checks that the run finishes, raises no traceback, and that the pipeline
reports the four results as stored.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest

from tests._support import PROJECT_ROOT

SCRIPT = "\n1\n2\n5\n1\n\n3\n\n0\n4\n0\nP\n\n0\n"


class MainMenuSmokeTest(unittest.TestCase):
    def test_orf_rna_protein_flow_is_stored(self):
        with tempfile.TemporaryDirectory() as log_dir:
            env = dict(os.environ, SYNAPSE_LOG_DIR=log_dir, PYTHONIOENCODING="utf-8")
            completed = subprocess.run(
                [sys.executable, "main.py"],
                cwd=PROJECT_ROOT,
                input=SCRIPT,
                capture_output=True,
                text=True,
                encoding="utf-8",
                timeout=120,
                env=env,
            )
        output = completed.stdout + completed.stderr
        self.assertEqual(completed.returncode, 0, output[-2000:])
        self.assertNotIn("Traceback", output)
        self.assertNotIn("ERROR", output.replace("ERRORS", ""))
        for label in ("Sequence Overview", "ORF Analysis", "RNA Analysis", "Protein Analysis"):
            self.assertRegex(output, rf"{label}\s+COMPLETED \(1 stored\)")
        self.assertIn("Selected ORF was used as RNA input", output)

    def test_report_export_from_menu(self):
        import json
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report.json"
            script = f"\n1\n2\n5\n1\n\nR\n{report}\n\n0\n"
            env = dict(os.environ, SYNAPSE_LOG_DIR=tmp, PYTHONIOENCODING="utf-8")
            completed = subprocess.run(
                [sys.executable, "main.py"], cwd=PROJECT_ROOT, input=script,
                capture_output=True, text=True, encoding="utf-8", timeout=120, env=env,
            )
            output = completed.stdout + completed.stderr
            self.assertEqual(completed.returncode, 0, output[-2000:])
            self.assertNotIn("Traceback", output)
            data = json.loads(report.read_text(encoding="utf-8"))
        modules = {entry["module"] for entry in data["analyses"]}
        self.assertIn("orf", modules)
        self.assertEqual(data["input"]["header"], "TEST_8_REFERENCE")
        orf_entry = next(e for e in data["analyses"] if e["module"] == "orf")
        self.assertEqual(orf_entry["parameters"]["min_length_nt"], 75)


if __name__ == "__main__":
    unittest.main()
