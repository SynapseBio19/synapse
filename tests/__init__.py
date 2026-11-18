"""Synapse automated test suite (standard-library ``unittest``).

Run from the project root:

    python -m unittest discover -s tests -t .

Tests never contact the network: ``tests._support`` blocks sockets.
"""

import os as _os

# The suite runs in a fixed environment: Synapse's per-user settings would
# otherwise change what it tests (found on 2026-11-19, when the owner set
# SYNAPSE_BLAST_DB_ROOT to a folder outside the project and the database-root
# checks failed). The original values stay available for the opt-in tests
# that need the user's real setup (tests/test_blast_refseq_rna.py).
USER_ENVIRONMENT = {}
for _name in ("SYNAPSE_BLAST_DB_ROOT", "SYNAPSE_HOME", "SYNAPSE_LOG_DIR", "SYNAPSE_NCBI_TOOL",
              "SYNAPSE_NCBI_EMAIL"):
    if _name in _os.environ:
        USER_ENVIRONMENT[_name] = _os.environ.pop(_name)

