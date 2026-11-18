"""Locating the NCBI BLAST+ programs and the database update script."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional

from analyses._blast.paths import ANALYSES_DIR
from analyses._blast.config import BLAST_EXECUTABLES, BLAST_UPDATE_SCRIPT
from analyses._blast.models import BlastAnalysisError


class BlastExecutableManager:
    """Locate BLAST+ executables without hard-coding a user-specific install path."""

    def __init__(self, blast_bin: Optional[str | Path] = None) -> None:
        self.explicit_bin = Path(blast_bin).expanduser() if blast_bin else None

    def _candidate_bins(self) -> List[Path]:
        candidates: List[Path] = []

        if self.explicit_bin:
            candidates.append(self.explicit_bin)

        for env_name in ("SYNAPSE_BLAST_BIN", "NCBI_BLAST_BIN"):
            value = os.environ.get(env_name, "").strip()
            if value:
                candidates.append(Path(value).expanduser())

        for env_name in ("BLAST_HOME", "NCBI_BLAST_HOME"):
            value = os.environ.get(env_name, "").strip()
            if value:
                root = Path(value).expanduser()
                candidates.extend([root / "bin", root])

        # Common Windows installation roots. These are only discovery hints;
        # no user-specific path is assumed.
        if os.name == "nt":
            for env_name in ("ProgramFiles", "ProgramFiles(x86)", "LOCALAPPDATA"):
                value = os.environ.get(env_name, "").strip()
                if not value:
                    continue
                base = Path(value)
                candidates.extend(
                    [
                        base / "NCBI" / "blast",
                        base / "NCBI",
                        base / "blast",
                    ]
                )
                try:
                    candidates.extend(
                        [path / "bin" for path in base.glob("NCBI/blast-*/bin")]
                    )
                except OSError:
                    pass

        # The executable may live next to this module in a portable Synapse
        # deployment or inside a local blast*/bin directory.
        module_root = ANALYSES_DIR
        candidates.extend(
            [
                module_root / "bin",
                module_root / "blast" / "bin",
            ]
        )
        try:
            candidates.extend([path / "bin" for path in module_root.glob("blast-*/bin")])
        except OSError:
            pass

        unique: List[Path] = []
        seen: set[str] = set()
        for candidate in candidates:
            key = str(candidate.resolve()) if candidate.exists() else str(candidate)
            if key not in seen:
                seen.add(key)
                unique.append(candidate)
        return unique

    def resolve(self, executable: str) -> Optional[str]:
        """Return a usable executable path, or None if unavailable."""
        direct = shutil.which(executable)
        if direct:
            return direct

        names = [executable]
        if os.name == "nt" and not executable.lower().endswith(".exe"):
            names.insert(0, executable + ".exe")

        for directory in self._candidate_bins():
            for name in names:
                candidate = directory / name
                if candidate.is_file():
                    return str(candidate)

        return None

    def available(self) -> Dict[str, Optional[str]]:
        return {name: self.resolve(name) for name in BLAST_EXECUTABLES}

    def require(self, executable: str) -> str:
        path = self.resolve(executable)
        if not path:
            raise BlastAnalysisError(
                f"BLAST+ executable '{executable}' was not found. "
                "Install NCBI BLAST+ or configure SYNAPSE_BLAST_BIN."
            )
        return path

    def resolve_perl(self) -> Optional[str]:
        return shutil.which("perl")

    def resolve_update_script(self) -> Optional[str]:
        """Locate NCBI's update_blastdb.pl without assuming a user-specific path.

        On Windows, ``blastn.exe`` may be discoverable through PATH while the
        accompanying Perl script is not itself on PATH.  In that common NCBI
        installation layout, the most reliable discovery method is to resolve
        a BLAST+ executable first and inspect its sibling ``bin`` directory.
        """
        direct = shutil.which(BLAST_UPDATE_SCRIPT)
        if direct:
            return direct

        # Prefer the directory containing a BLAST+ executable already visible
        # to the operating system.  This handles standard Windows installs
        # such as ``...\\NCBI\\blast-2.17.0+\\bin`` even when the Perl
        # script itself is not included in PATH.
        for executable in BLAST_EXECUTABLES:
            resolved_executable = shutil.which(executable)
            if not resolved_executable:
                continue
            candidate = Path(resolved_executable).resolve().parent / BLAST_UPDATE_SCRIPT
            if candidate.is_file():
                return str(candidate)

        for directory in self._candidate_bins():
            candidate = directory / BLAST_UPDATE_SCRIPT
            if candidate.is_file():
                return str(candidate)

        return None
