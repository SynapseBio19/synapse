"""The feature inventory names only code that exists (Phase 7, 2026-11-17).

`docs/FEATURE_STATUS.md` lists, for each feature, the functions that
implement it. On 2026-11-17 two of those names had been deleted months
earlier (D16 part 2) and several rows described a state the code had left
behind. This test keeps the "Function" column honest: every identifier in
it must occur in the application source, and every file it names must
exist. Execution statuses are checked against line coverage by hand (see
section 18 of the document); that needs a coverage run and is not part of
the suite.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCUMENT = ROOT / "docs" / "FEATURE_STATUS.md"
# analyses/**: the internal packages of Phase 7.1 (analyses/_blast, analyses/_primer).
SOURCES = [*ROOT.glob("analyses/**/*.py"), *ROOT.glob("core/*.py"), *ROOT.glob("ui/*.py"),
           *ROOT.glob("fileio/*.py"), *ROOT.glob("data/*.py"), ROOT / "main.py"]
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def function_cells():
    """(section, feature, cell) for every row of sections 1-15 that has a function column."""
    section, header, cells = None, None, []
    for line in DOCUMENT.read_text(encoding="utf-8").splitlines():
        heading = re.match(r"^## (\d+)\. ", line)
        if heading:
            section, header = int(heading.group(1)), None
            continue
        if section is None or section > 15 or not line.startswith("|"):
            continue
        row = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if set("".join(row)) <= set("-: "):
            continue
        if header is None:
            header = row
            continue
        values = dict(zip(header, row))
        cell = values.get("Function") or values.get("Function/class")
        if cell:
            cells.append((section, values.get("Feature") or values.get("Component") or "", cell))
    return cells


class FeatureInventoryNamesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = "\n".join(path.read_text(encoding="utf-8") for path in SOURCES)
        cls.identifiers = set(IDENTIFIER.findall(cls.source))

    def _exists(self, token):
        """A file path must exist; an identifier must occur in the source."""
        bare = token.strip().strip("()")
        if bare.endswith((".py", ".json")) or bare in {"data", "ui", "blast_results", "tests"}:
            candidates = [ROOT / bare, ROOT / "data" / bare, ROOT / "ui" / bare, ROOT / "blast_results" / bare]
            return any(path.exists() for path in candidates)
        if bare.startswith(("analyses.", "core.", "ui.", "fileio.", "data.")) and "(" not in bare:
            parts = bare.split(".")
            module = ROOT.joinpath(*parts[:2]).with_suffix(".py")
            if len(parts) == 2:
                return module.exists()
            bare = ".".join(parts[2:])
        wildcard = "*" in bare or "…" in bare or re.fullmatch(r"_v\d+", bare) is not None
        for name in IDENTIFIER.findall(bare):
            if len(name) < 3 or name in {"py", "main"}:
                continue
            found = name in self.source if wildcard else name in self.identifiers
            if not found:
                return False
        return True

    def test_every_named_function_exists(self):
        missing = []
        for section, feature, cell in function_cells():
            for token in re.findall(r"`([^`]+)`", cell):
                for part in re.split(r",\s*", token):
                    if part and not self._exists(part):
                        missing.append(f"§{section} {feature}: {part}")
        self.assertEqual(missing, [], "FEATURE_STATUS names code that does not exist")

    def test_the_check_detects_a_removed_name(self):
        self.assertFalse(self._exists("_effective_monovalent_salt_mM"))
        self.assertFalse(self._exists("_professional_failure_diagnostics"))
        self.assertTrue(self._exists("calculate_salt_correction"))


if __name__ == "__main__":
    unittest.main()
