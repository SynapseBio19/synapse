"""Fail on any pyright error except the documented ones (Phase 7, 2026-11-17).

Runs `pyright --pythonversion 3.10 --outputjson` with the current
interpreter and compares the errors with an allow-list. The allowed errors
are the findings in `analyses/crispr.py`, which Phase 7 does not edit (owner
decision of 2026-11-17: CRISPR work is excluded). Any other error, or more
errors in that file than recorded, fails the check, so the count can only
go down. Usage: `python scripts/ci/pyright_gate.py`.
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter
from pathlib import PurePath

# Errors that are known and accepted, by file, with the reason.
ALLOWED = {
    "analyses/crispr.py": 11,   # CRISPR module excluded from changes in Phase 7
}


def main() -> int:
    command = [sys.executable, "-m", "pyright", "--pythonversion", "3.10",
               "--pythonpath", sys.executable, "--outputjson"]
    completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8")
    try:
        report = json.loads(completed.stdout)
    except json.JSONDecodeError:
        print(completed.stdout or completed.stderr)
        print("ERROR: pyright did not produce a JSON report.")
        return 2
    errors = [item for item in report.get("generalDiagnostics", []) if item.get("severity") == "error"]
    by_file: Counter[str] = Counter()
    for item in errors:
        parts = PurePath(item["file"]).parts
        anchor = max((i for i, part in enumerate(parts) if part in {"analyses", "core", "ui", "fileio",
                                                                      "data", "tests"}), default=None)
        relative = "/".join(parts[anchor:]) if anchor is not None else PurePath(item["file"]).name
        by_file[relative] += 1
        if relative not in ALLOWED:
            line = item["range"]["start"]["line"] + 1
            print(f"{relative}:{line}: {item.get('rule', '')}: {item['message'].splitlines()[0]}")
    failures = []
    for relative, count in sorted(by_file.items()):
        allowed = ALLOWED.get(relative, 0)
        if count > allowed:
            failures.append(f"{relative}: {count} error(s), {allowed} allowed")
    print(f"pyright: {len(errors)} error(s); allowed by file: {ALLOWED}")
    for relative, allowed in ALLOWED.items():
        if by_file.get(relative, 0) < allowed:
            print(f"NOTE: {relative} has {by_file.get(relative, 0)} error(s), fewer than the {allowed} "
                  "allowed; lower the number in scripts/ci/pyright_gate.py.")
    if failures:
        print("ERROR: new pyright errors:\n  " + "\n  ".join(failures))
        return 1
    print("OK: no pyright error outside the documented allow-list.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
