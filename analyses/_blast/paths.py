"""Locations Synapse resolves relative to its installation."""

from __future__ import annotations

import os
from pathlib import Path

# The ``analyses`` directory. Before Phase 7.1 the BLAST code lived in
# analyses/blast.py and used ``Path(__file__).resolve().parent`` for it.
ANALYSES_DIR = Path(__file__).resolve().parent.parent


def _resolve_synapse_user_path(path: str | Path) -> Path:
    """Resolve a user-supplied path independently of the process CWD.

    Synapse can be launched from Windows Terminal, VS Code, a shortcut, or
    another application. In some launch modes Windows uses ``C:\\Windows\\System32``
    as the current working directory. Relative paths entered by the user
    should nevertheless be interpreted relative to the Synapse project root.

    Absolute paths are preserved. Relative paths are first interpreted from
    the Synapse project root (the directory containing ``analyses``).
    """
    candidate = Path(path).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()

    working_directory = Path.cwd()
    # Windows shortcuts can launch console applications with System32 as the
    # process CWD. In that special case, use Synapse's application root as a
    # practical fallback; otherwise preserve normal relative-path semantics.
    if os.name == "nt" and working_directory.resolve() == Path(os.environ.get("SystemRoot", r"C:\Windows")).resolve() / "System32":
        working_directory = ANALYSES_DIR.parent
    return (working_directory / candidate).resolve()
