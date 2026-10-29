"""Golden-master comparison for every analysis module.

Fails when the output of a module changes. To accept an intended change:

    SYNAPSE_UPDATE_GOLDEN=1 python -m unittest tests.test_golden_master

and record the reason in docs/CHANGES.md.
"""

from __future__ import annotations

import json
import unittest

from tests import _support
from tests.golden_scenarios import SCENARIOS


def _make_test(name, builder):
    def test(self):
        data = _support.normalize_for_json(builder())
        # JSON round-trip turns tuples into lists, as in the stored file.
        data = json.loads(json.dumps(data))
        if _support.update_golden_requested() or not _support.golden_path(name).exists():
            _support.write_golden(name, data)
            self.skipTest(f"golden snapshot '{name}' written")
        expected = _support.load_golden(name)
        self.maxDiff = 4000
        self.assertEqual(expected, data, f"Golden master '{name}' changed")
    test.__name__ = f"test_golden_{name}"
    return test


class GoldenMasterTests(unittest.TestCase):
    pass


for _name, _builder in SCENARIOS.items():
    setattr(GoldenMasterTests, f"test_golden_{_name}", _make_test(_name, _builder))


if __name__ == "__main__":
    unittest.main()
