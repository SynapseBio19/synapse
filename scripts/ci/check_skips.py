"""Summarize skipped tests from a `unittest -v` log and enforce what may skip.

`python scripts/ci/check_skips.py LOG` lists every skipped test with its
reason (also into the GitHub step summary when `GITHUB_STEP_SUMMARY` is
set) and always succeeds.

`python scripts/ci/check_skips.py LOG --only "set SYNAPSE_BENCHMARK=1"`
fails if any test skipped for a reason that does not contain one of the
given texts. The reference job uses it so that a missing reference package
cannot pass unnoticed.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import Counter

# A test is shown by its id or, if it has a docstring, by the docstring's first
# line; unittest quotes the reason with double quotes when it contains an apostrophe.
SKIP = re.compile(r"^(?P<test>.*?) \.\.\. skipped (?P<quote>['\"])(?P<reason>.*)(?P=quote)$")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("log")
    parser.add_argument("--only", action="append", default=[],
                        help="reason text that is allowed to skip (repeatable)")
    args = parser.parse_args()
    with open(args.log, encoding="utf-8", errors="replace") as handle:
        skips = [m.groupdict() for m in (SKIP.match(line.rstrip("\n")) for line in handle) if m]
    reasons = Counter(skip["reason"] for skip in skips)
    lines = [f"{len(skips)} skipped test(s)"] + [f"- {count} x {reason}" for reason, count in reasons.most_common()]
    print("\n".join(lines))
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as handle:
            handle.write("### Skipped tests\n\n" + "\n".join(lines) + "\n\n")
    if args.only:
        unexpected = [skip for skip in skips if not any(text in skip["reason"] for text in args.only)]
        if unexpected:
            print("ERROR: tests skipped for a reason that this job does not allow:")
            for skip in unexpected:
                print(f"  {skip['test']}: {skip['reason']}")
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
