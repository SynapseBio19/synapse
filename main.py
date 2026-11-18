"""Start Synapse from the source folder: ``python main.py``.

The session itself is ``synapse_gen/app.py`` (moved there in Phase 9 so that
it can be installed; the installed command is ``synapse``). Importing
``main`` returns that module itself, so ``main.main()``, ``main.current_folder``
and patches of ``main`` keep working for scripts and tests.
"""

import sys

from synapse_gen import app as _app
from synapse_gen.app import *  # noqa: F401,F403  (the names, for static analysis)

if __name__ == "__main__":
    raise SystemExit(_app.cli())
else:
    sys.modules[__name__] = _app
