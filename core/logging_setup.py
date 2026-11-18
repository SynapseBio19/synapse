"""Central logging configuration for the interactive application.

Log file: ``logs/synapse.log`` in Synapse's user folder: the project folder
when Synapse runs from source, a per-user folder when it is installed
(``core.paths``). Override with the environment variable ``SYNAPSE_LOG_DIR``. Analysis modules log through
``logging.getLogger("synapse.<module>")``; they never configure handlers
themselves, so library use and tests stay silent.
"""

from __future__ import annotations

import logging
import os
from logging.handlers import RotatingFileHandler
from pathlib import Path

from core.paths import SOURCE_ROOT, user_data_dir

PROJECT_ROOT = SOURCE_ROOT


def log_directory() -> Path:
    configured = os.environ.get("SYNAPSE_LOG_DIR", "").strip()
    return Path(configured).expanduser() if configured else user_data_dir() / "logs"


def configure_logging(level: int = logging.INFO) -> Path:
    """Attach a rotating file handler to the 'synapse' logger; return the log path."""
    directory = log_directory()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "synapse.log"

    logger = logging.getLogger("synapse")
    logger.setLevel(level)
    if not any(isinstance(h, RotatingFileHandler) for h in logger.handlers):
        handler = RotatingFileHandler(path, maxBytes=1_000_000, backupCount=3, encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
        logger.addHandler(handler)
    logger.propagate = False
    return path
