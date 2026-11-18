"""Source rules for messages and prompts (Phase 7, package G).

Phase 7 found five hand-padded "Label:     value" lines in the single-primer
analysis although TERMINAL_UI_STANDARD said none remained, and three BLAST
prompts that began with a carriage return. These rules keep both from
returning. The CRISPR dialog is outside Phase 7 and is not checked here.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FILES = [path for path in [*ROOT.glob("ui/*.py"), ROOT / "synapse_gen" / "app.py"] if path.name != "crispr_menu.py"]

HAND_PADDED = re.compile(r"""print\(f?["'][^"'\n]{1,40}:\s{3,}""")
PROMPT_WITH_BREAK = re.compile(r"""\bask\(\s*f?["']\\[rn]""")


def findings(pattern):
    hits = []
    for path in FILES:
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            if pattern.search(line):
                hits.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()[:80]}")
    return hits


class MessageRuleTests(unittest.TestCase):
    def test_fields_are_drawn_with_style_field(self):
        self.assertEqual(findings(HAND_PADDED), [])

    def test_prompts_do_not_begin_with_a_line_break(self):
        self.assertEqual(findings(PROMPT_WITH_BREAK), [])

    def test_the_rules_detect_the_old_forms(self):
        self.assertTrue(HAND_PADDED.search('print(f"ΔH:                     {x:.2f} kcal/mol")'))
        self.assertTrue(PROMPT_WITH_BREAK.search("name = style.ask('\\r\\nDatabase name (0 cancel)')"))
        self.assertFalse(HAND_PADDED.search("style.field('ΔH', value)"))


if __name__ == "__main__":
    unittest.main()
