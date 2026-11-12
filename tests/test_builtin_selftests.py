"""Run every built-in module self-test from the automated suite."""

from __future__ import annotations

import contextlib
import importlib
import io
import unittest

from tests import _support  # noqa: F401  (network + input guard)

SELF_TESTS = [
    ("analyses.sequence", "run_sequence_self_test"),
    ("analyses.restriction", "run_restriction_self_test"),
    ("analyses.cloning", "run_cloning_self_test"),
    ("ui.cloning_menu", "run_cloning_ui_self_test"),
    ("analyses.mutation", "run_mutation_self_test"),
    ("analyses.variant_analysis", "run_variant_analysis_self_test"),
    ("analyses.primer", "run_primer_self_test"),
    ("analyses.primer", "run_thermodynamic_self_test"),
    ("analyses.primer", "run_primer_v3_self_test_final"),
    ("analyses.primer", "run_primer_v4_screening_self_test"),
    ("analyses.crispr", "run_self_test"),
    ("analyses.blast", "run_self_test"),
    ("analyses.literature", "self_test"),
]


def _run(module_name, function_name):
    module = importlib.import_module(module_name)
    function = getattr(module, function_name)
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            return function(verbose=False)
        except TypeError:
            return function()


class BuiltinSelfTests(unittest.TestCase):
    def test_all_builtin_self_tests_pass(self):
        for module_name, function_name in SELF_TESTS:
            with self.subTest(f"{module_name}.{function_name}"):
                result = _run(module_name, function_name)
                if isinstance(result, dict):
                    failed = [k for k, v in result.items() if v is not True]
                    self.assertEqual([], failed, f"failing checks: {failed}")
                else:
                    self.assertIs(result, True)


if __name__ == "__main__":
    unittest.main()
