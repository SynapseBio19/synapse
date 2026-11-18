"""Runtime views of local databases: command references, temporary alias
files for nested volumes, and the BLAST process environment."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, List, Tuple

from analyses._blast.models import BlastAnalysisError, LocalBlastDatabase


def _build_local_database_command_reference(database_path: str | Path) -> Tuple[str, Path]:
    """Return a BLAST database reference and working directory.

    BLAST+ on Windows is generally path-safe when launched through
    ``subprocess`` with an argument list, but some BLAST+ builds and database
    layouts behave more reliably when the database prefix is supplied relative
    to its own directory. This also avoids ambiguity for paths containing
    spaces (for example ``Project Synapse``).

    The returned reference is only the database prefix name; the caller must
    set ``cwd`` to the returned parent directory. No shell is used.
    """
    resolved = Path(database_path).expanduser().resolve()
    parent = resolved.parent
    return resolved.name, parent


def _volume_prefix_from_stored_path(volume_path: str | Path) -> Path:
    """Resolve a stored volume path to the actual BLAST database prefix.

    Nested NCBI extraction layouts are stored by Synapse as directories for
    display purposes, e.g. ``.../refseq_rna.01``. The BLAST prefix inside that
    directory is ``.../refseq_rna.01/refseq_rna.01``. Flat prefixes are kept as
    they are.
    """
    resolved = Path(volume_path).expanduser().resolve()
    if resolved.is_dir():
        nested_prefix = resolved / resolved.name
        extensions = (".nhr", ".nin", ".nsq", ".phr", ".pin", ".psq")
        if any(Path(str(nested_prefix) + extension).is_file() for extension in extensions):
            return nested_prefix
    return resolved


def _create_runtime_blast_alias(
    database: LocalBlastDatabase,
    runtime_directory: str | Path,
) -> Tuple[str, Path]:
    """Create a temporary flat BLAST view for a nested multi-volume database.

    NCBI's preformatted multi-volume databases are normally exposed to BLAST
    as a flat set of files:

        refseq_rna.nal
        refseq_rna.00.nhr / .nin / .nsq / ...
        refseq_rna.01.nhr / .nin / .nsq / ...
        ...

    Synapse also supports the extraction layout produced by manually unpacking
    each NCBI archive into a volume directory, for example::

        refseq_rna/refseq_rna.00/refseq_rna.00.nhr
        refseq_rna/refseq_rna.01/refseq_rna.01.nhr
        ...

    BLAST+ does not treat that nested layout as the canonical multi-volume
    database when it is presented through a synthetic alias. In particular,
    v5 metadata such as ``refseq_rna.ndb`` is then interpreted incorrectly for
    later volumes. The robust solution is therefore a temporary *flat view*.

    The view uses hard links, not file copies. Hard links consume essentially
    no additional database storage and leave the original NCBI files untouched.
    The runtime directory is on the same filesystem as the source database so
    Windows can create the hard links. The alias keeps the canonical database
    name and volume names, which lets BLAST+ use its normal multi-volume logic.
    """
    requested_runtime = Path(runtime_directory).expanduser().resolve()
    database_path = Path(database.path).expanduser().resolve()
    requested_runtime.mkdir(parents=True, exist_ok=True)

    extension = ".pal" if database.database_type == "protein" else ".nal"
    logical_name = Path(database.name).name
    logical_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", logical_name) or "synapse_database"

    # A BLAST multi-volume alias must use the same logical root as its volume
    # prefixes (for example refseq_rna + refseq_rna.00 ... refseq_rna.20).
    # This is intentionally NOT a *_runtime.nal alias.
    alias_path = requested_runtime / f"{logical_name}{extension}"

    volume_prefixes = [
        _volume_prefix_from_stored_path(path)
        for path in database.volume_paths
    ]
    if not volume_prefixes:
        raise BlastAnalysisError(
            f"Cannot create runtime BLAST view for '{database.name}': "
            "no volume paths were detected."
        )

    # The source database must be on the same filesystem because Windows
    # hard-links cannot cross volumes. tempfile.TemporaryDirectory callers
    # already place the directory next to the database for this reason.
    try:
        requested_runtime.relative_to(database_path.parent.parent)
    except ValueError:
        # This is only a diagnostic guard; a same-volume check is performed by
        # os.link below and produces the actionable error if needed.
        pass

    db_entries: List[str] = []
    created_links: List[Path] = []

    try:
        for volume_prefix in volume_prefixes:
            volume_name = volume_prefix.name
            if not re.fullmatch(rf"{re.escape(logical_name)}\.\d+", volume_name, re.IGNORECASE):
                raise BlastAnalysisError(
                    f"Unexpected BLAST volume name '{volume_name}' for database '{logical_name}'."
                )

            # Reproduce the canonical flat BLAST volume prefix. Link every
            # source file belonging to this volume, including v5 metadata on
            # the first volume when present.
            source_parent = volume_prefix.parent
            pattern = f"{volume_name}.*"
            source_files = sorted(source_parent.glob(pattern))
            if not source_files:
                raise BlastAnalysisError(
                    f"No files found for BLAST volume '{volume_prefix}'."
                )

            for source_file in source_files:
                if not source_file.is_file():
                    continue
                target = requested_runtime / source_file.name
                try:
                    if target.exists() or target.is_symlink():
                        target.unlink()
                    os.link(str(source_file), str(target))
                    created_links.append(target)
                except OSError as exc:
                    raise BlastAnalysisError(
                        "Could not create a temporary hard link for BLAST database "
                        f"file '{source_file}': {exc}"
                    ) from exc

            db_entries.append(volume_name)

        # The original NCBI alias contains the database-level v5 files in the
        # first volume directory. Those files do not carry a volume suffix and
        # therefore are not covered by the volume glob above.
        first_volume_dir = volume_prefixes[0].parent
        for source_file in sorted(first_volume_dir.iterdir()):
            if not source_file.is_file():
                continue
            if source_file.name == alias_path.name:
                continue
            # Only link database-level metadata that BLAST may need. Taxonomy
            # files are included as well so the runtime view remains complete.
            if re.fullmatch(rf"{re.escape(logical_name)}\.\d+\..+", source_file.name, re.IGNORECASE):
                continue
            if source_file.name.lower() in {
                "taxdb.btd",
                "taxdb.bti",
                "taxonomy4blast.sqlite3",
            } or source_file.name.lower() in {
                f"{logical_name}.ndb",
                f"{logical_name}.nos",
                f"{logical_name}.not",
                f"{logical_name}.ntf",
                f"{logical_name}.nto",
            }:
                target = requested_runtime / source_file.name
                try:
                    if target.exists() or target.is_symlink():
                        target.unlink()
                    os.link(str(source_file), str(target))
                    created_links.append(target)
                except OSError as exc:
                    raise BlastAnalysisError(
                        f"Could not create a temporary hard link for '{source_file}': {exc}"
                    ) from exc

        title = database.title or f"Synapse runtime view for {database.name}"
        alias_path.write_text(
            "# Synapse temporary BLAST alias\r\n"
            f"TITLE {title}\r\n"
            f"DBLIST {' '.join(db_entries)}\r\n",
            encoding="utf-8",
            newline="",
        )

        return alias_path.stem, alias_path.parent

    except Exception:
        # Never leave a partially constructed runtime view behind.
        for path in reversed(created_links):
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass
        try:
            alias_path.unlink(missing_ok=True)
        except OSError:
            pass
        raise

def _needs_runtime_blast_alias(database: LocalBlastDatabase) -> bool:
    """Return whether the physical database layout needs a runtime alias."""
    if database.volume_count <= 1 or not database.volume_paths:
        return False

    database_path = Path(database.path).expanduser().resolve()
    parent_name = database_path.parent.name
    return bool(re.fullmatch(r".+\.\d+", parent_name))


def _build_local_blast_environment(database_path: str | Path) -> Dict[str, str]:
    """Build a BLAST+ environment with an explicit database search path."""
    env = os.environ.copy()
    resolved = Path(database_path).expanduser().resolve()
    parent = resolved if resolved.is_dir() else resolved.parent
    existing = env.get("BLASTDB", "").strip()
    parent_text = str(parent)
    if existing:
        env["BLASTDB"] = os.pathsep.join([parent_text, existing])
    else:
        env["BLASTDB"] = parent_text
    return env
