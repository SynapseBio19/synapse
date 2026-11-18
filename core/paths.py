"""Where Synapse reads and writes its files (Phase 9, 2026-11-19).

Two ways of running Synapse:

- **From the source folder** (``python main.py``, or an editable install):
  ``pyproject.toml`` sits next to the packages. Logs, exports and the
  default BLAST database root stay in the project folder, as before.
- **Installed** (``pip install synapse-gen``, command ``synapse``): the
  packages live in ``site-packages``, which is not a place for user files.
  Synapse then uses a per-user folder: ``%LOCALAPPDATA%\\Synapse`` on
  Windows, ``$XDG_DATA_HOME/synapse`` (``~/.local/share/synapse``)
  elsewhere.

``SYNAPSE_HOME`` overrides both. The more specific variables
(``SYNAPSE_LOG_DIR``, ``SYNAPSE_BLAST_DB_ROOT``) still take precedence for
their own files.
"""

from __future__ import annotations

import os
from pathlib import Path

HOME_ENV = "SYNAPSE_HOME"
SOURCE_ROOT = Path(__file__).resolve().parent.parent


def is_source_checkout() -> bool:
    """True when Synapse runs from its project folder, not from site-packages."""
    return (SOURCE_ROOT / "pyproject.toml").is_file()


def user_data_dir() -> Path:
    """The folder for logs, exports and (by default) BLAST databases."""
    configured = os.environ.get(HOME_ENV, "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    if is_source_checkout():
        return SOURCE_ROOT
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA", "").strip()
        return (Path(base) if base else Path.home() / "AppData" / "Local") / "Synapse"
    base = os.environ.get("XDG_DATA_HOME", "").strip()
    return (Path(base) if base else Path.home() / ".local" / "share") / "synapse"


def example_sequence_file() -> Path:
    """The bundled example FASTA (installed with the ``data`` package)."""
    import data

    return Path(data.__file__).resolve().parent / "examples" / "dna.fasta"
