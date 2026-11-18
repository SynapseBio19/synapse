"""The local BLAST database manager: discovery, verification, manifest,
NCBI download, custom databases and removal."""

from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from analyses._blast.config import (
    DEFAULT_LOCAL_DB_ROOT,
    LOCAL_DATABASE_CATALOG,
    LOCAL_DB_MANIFEST_NAME,
    LOCAL_DB_ROOT_ENV,
)
from analyses._blast.models import BlastAnalysisError, LocalBlastDatabase, _utc_now
from analyses._blast.executables import BlastExecutableManager
from analyses._blast.runtime import (
    _build_local_blast_environment,
    _build_local_database_command_reference,
    _create_runtime_blast_alias,
    _needs_runtime_blast_alias,
)


def _parse_blastdbcmd_info(output: str) -> Dict[str, Any]:
    """Extract stable summary fields from blastdbcmd -info output."""
    text = output or ""
    info: Dict[str, Any] = {
        "title": None,
        "sequence_count": None,
        "total_length": None,
        "blastdb_version": None,
        "database_date": None,
        "volumes": [],
    }

    date_match = re.search(r"^\s*Date:\s*(.+?)(?:\t|\s{2,}Longest|$)", text, re.MULTILINE)
    if date_match:
        info["database_date"] = date_match.group(1).strip()

    database_match = re.search(r"^Database:\s*(.+)$", text, re.MULTILINE)
    if database_match:
        info["title"] = database_match.group(1).strip()

    count_match = re.search(
        r"([\d,]+)\s+sequences;\s+([\d,]+)\s+total\s+bases",
        text,
        re.IGNORECASE,
    )
    if count_match:
        info["sequence_count"] = int(count_match.group(1).replace(",", ""))
        info["total_length"] = int(count_match.group(2).replace(",", ""))

    version_match = re.search(
        r"BLASTDB\s+Version:\s*([0-9]+)",
        text,
        re.IGNORECASE,
    )
    if version_match:
        info["blastdb_version"] = version_match.group(1)

    in_volumes = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.lower() == "volumes:":
            in_volumes = True
            continue
        if in_volumes and stripped:
            if stripped.startswith("Database:") or stripped.startswith("Date:"):
                break
            info["volumes"].append(stripped)

    return info


def _extract_blast_volume_index(volume_name: str) -> Optional[int]:
    """Return the numeric BLAST volume index while preserving zero padding elsewhere.

    NCBI alias files and extracted database layouts can represent the same
    volume as ``refseq_rna.0`` or ``refseq_rna.00``. The numeric index is the
    stable identity; the actual filesystem name (including zero padding) is
    preserved by the resolver.
    """
    match = re.search(r"\.(\d+)$", str(volume_name).strip())
    if not match:
        return None
    try:
        return int(match.group(1))
    except ValueError:
        return None


def _parse_blast_alias_volume_names(alias_path: Path) -> List[str]:
    """Read volume names declared by an NCBI BLAST alias (.nal/.pal) file.

    NCBI multi-volume databases can be stored in dedicated folders such as
    ``refseq_rna.00/`` while the alias file lives one directory above them.
    The alias remains the authoritative source for which volumes belong to the
    logical database, so Synapse uses DBLIST when it is available.
    """
    if not alias_path.is_file():
        return []

    try:
        payload = alias_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []

    match = re.search(
        r"(?ims)^\s*DBLIST\s+(.+?)(?=^\s*[A-Z][A-Z0-9_]*\s+|\Z)",
        payload,
    )
    if not match:
        return []

    block = match.group(1)
    values = re.findall(r'"([^"]+)"', block)
    if not values:
        values = re.findall(r"[^\s]+", block)

    normalized: List[str] = []
    seen: set[str] = set()
    for value in values:
        value = value.strip().strip('"').replace("\\", "/")
        if not value:
            continue
        name = Path(value).name
        if name and name not in seen:
            seen.add(name)
            normalized.append(name)
    return normalized


class LocalBlastDatabaseManager:
    """
    Manage NCBI preformatted BLAST databases and custom FASTA databases.

    The manager never downloads a database implicitly. NCBI's
    update_blastdb.pl is used for official preformatted database downloads,
    while makeblastdb is used for custom FASTA collections.
    """

    def __init__(
        self,
        root: Optional[str | Path] = None,
        executable_manager: Optional[BlastExecutableManager] = None,
    ) -> None:
        configured = os.environ.get(LOCAL_DB_ROOT_ENV, "").strip()
        self.root = Path(root or configured or DEFAULT_LOCAL_DB_ROOT).expanduser()
        self.root = self.root.resolve()
        self.executable_manager = executable_manager or BlastExecutableManager()
        self.manifest_path = self.root / LOCAL_DB_MANIFEST_NAME

    def ensure_root(self) -> Path:
        self.root.mkdir(parents=True, exist_ok=True)
        return self.root

    def _load_manifest(self) -> Dict[str, Any]:
        if not self.manifest_path.is_file():
            return {}
        try:
            payload = json.loads(self.manifest_path.read_text(encoding="utf-8"))
            return payload if isinstance(payload, dict) else {}
        except (OSError, json.JSONDecodeError):
            return {}

    def _save_manifest(self, payload: Dict[str, Any]) -> None:
        self.ensure_root()
        temp_path = self.manifest_path.with_suffix(".tmp")
        temp_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        temp_path.replace(self.manifest_path)

    def _safe_name(self, name: str) -> str:
        """Sanitize a database name while allowing safe subdirectories."""
        raw = str(name).strip().replace("\\", "/")
        if not raw:
            raise BlastAnalysisError("Database name is empty.")

        parts: List[str] = []
        for part in raw.split("/"):
            part = part.strip()
            if not part or part in {".", ".."}:
                if part == "..":
                    raise BlastAnalysisError("Database path traversal is not allowed.")
                continue
            cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", part)
            cleaned = cleaned.strip("._")
            if cleaned in {"", ".", ".."}:
                raise BlastAnalysisError("Invalid database name.")
            parts.append(cleaned)

        if not parts:
            raise BlastAnalysisError("Database name is empty after sanitization.")
        return "/".join(parts)

    def _prefix(self, name: str) -> Path:
        """Return a database prefix guaranteed to remain below self.root."""
        safe_name = self._safe_name(name)
        prefix = (self.root / Path(*safe_name.split("/"))).resolve()
        try:
            prefix.relative_to(self.root)
        except ValueError as exc:
            raise BlastAnalysisError("Database path escapes the configured BLAST root.") from exc
        return prefix

    def _expected_extensions(self, database_type: str) -> Tuple[str, ...]:
        # .njs/.pjs: the JSON metadata file of BLAST+ version 5 databases
        # (Phase 7: without it, Remove Database left that file behind).
        if database_type == "protein":
            return (".pal", ".phr", ".pin", ".psq", ".pdb", ".pog", ".pos", ".pot", ".ptf", ".pto", ".pni",
                    ".pjs")
        return (".nal", ".nhr", ".nin", ".nsq", ".ndb", ".nni", ".nog", ".nos", ".not", ".ntf", ".nto",
                ".njs")

    def _find_existing_prefix(self, name: str, database_type: str) -> Path:
        """Resolve a database prefix in flexible local BLAST layouts.

        Supported layouts include:

        1. Flat databases::

               blast_db/testdb.nhr
               blast_db/testdb.nin
               blast_db/testdb.nsq

        2. Database folders::

               blast_db/refseq_rna/refseq_rna.nal
               blast_db/refseq_rna/refseq_rna.00/
               blast_db/refseq_rna/refseq_rna.01/

        3. Nested volume files::

               blast_db/refseq_rna/refseq_rna.00.nhr

        The logical database name remains the requested name; the returned
        path is the physical BLAST prefix used by BLAST+.
        """
        direct = self._prefix(name)
        extensions = self._expected_extensions(database_type)
        if any(Path(str(direct) + extension).is_file() for extension in extensions):
            return direct

        # A database container may have the same name as the logical database,
        # with the actual alias prefix one level below it.
        if direct.is_dir():
            nested_prefix = direct / direct.name
            if any(
                Path(str(nested_prefix) + extension).is_file()
                for extension in extensions
            ):
                return nested_prefix

            # NCBI archives can be extracted volume-by-volume. In that case
            # the first archive may contain the logical alias and shared BLAST
            # metadata inside the ``.00`` volume directory, while ``.01`` ...
            # ``.20`` remain sibling volume directories. The alias is still
            # the logical database; it must not be rejected merely because its
            # physical files live in the first volume directory.
            volume_dir_pattern = re.compile(r"^.+\.\d+$", re.IGNORECASE)
            try:
                alias_candidate = direct / f"{direct.name}.nal"
                protein_alias_candidate = direct / f"{direct.name}.pal"
                if alias_candidate.is_file() or protein_alias_candidate.is_file():
                    return direct

                first_volume = direct / f"{direct.name}.00"
                if first_volume.is_dir():
                    for extension in extensions:
                        candidate = first_volume / f"{direct.name}{extension}"
                        if candidate.is_file():
                            return first_volume / direct.name

                    # More robustly locate the alias in the first numbered
                    # volume directory without treating it as a separate DB.
                    for marker in first_volume.iterdir():
                        if not marker.is_file():
                            continue
                        if marker.name.lower() in {
                            f"{direct.name}.nal".lower(),
                            f"{direct.name}.pal".lower(),
                        }:
                            return first_volume / direct.name
            except OSError:
                pass

        # Search recursively for the requested database prefix. This supports
        # user-chosen roots and nested database containers without hard-coding
        # a Project Synapse path.
        leaf = direct.name
        try:
            matches: List[Path] = []
            volume_dir_pattern = re.compile(r"^.+\.\d+$", re.IGNORECASE)
            for extension in extensions:
                for path in self.root.rglob(f"{leaf}{extension}"):
                    if not path.is_file():
                        continue
                    # Never resolve a logical database prefix from inside a
                    # numbered BLAST volume directory. Those files belong to
                    # the parent database and are not a separate database.
                    if any(volume_dir_pattern.fullmatch(part) for part in path.parent.parts):
                        continue
                    matches.append(path.with_suffix(""))
            if matches:
                matches = sorted(
                    {path.resolve() for path in matches},
                    key=lambda path: (
                        0 if path.parent.name == leaf else 1,
                        len(path.parts),
                        str(path).lower(),
                    ),
                )
                return matches[0]
        except OSError:
            pass

        return direct

    def _volume_container_prefix(self, prefix: Path) -> Path:
        """Return the logical container prefix used to discover sibling volumes.

        Some manually extracted NCBI database archives place the logical alias
        and shared metadata inside ``name.00/`` while the remaining volumes are
        siblings of that directory. In that layout the physical alias prefix is
        ``name.00/name``, but volume discovery must happen relative to
        ``parent/name``.
        """
        base = prefix.name
        parent_name = prefix.parent.name
        if re.fullmatch(rf"{re.escape(base)}\.\d+", parent_name, re.IGNORECASE):
            return prefix.parent.parent / base
        return prefix

    def _discover_volume_prefixes(
        self,
        prefix: Path,
        database_type: str,
    ) -> Dict[str, Path]:
        """Discover BLAST volume prefixes below the alias directory.

        This deliberately handles both ``name.00.nhr`` files and the folder
        layout ``name.00/name.00.nhr`` used by extracted NCBI database archives.
        """
        extensions = self._expected_extensions(database_type)
        base = prefix.name
        volume_pattern = re.compile(rf"^{re.escape(base)}\.(\d+)$", re.IGNORECASE)
        found: Dict[str, Path] = {}

        def add_candidate(candidate: Path) -> None:
            match = volume_pattern.fullmatch(candidate.name)
            if match:
                found[candidate.name] = candidate.resolve()

        container = prefix.parent

        # Flat multi-volume files: name.00.nhr / name.01.nhr / ...
        try:
            for path in container.glob(f"{base}.*"):
                if path.is_file():
                    for extension in extensions:
                        if path.name.lower().endswith(extension.lower()):
                            add_candidate(path.with_suffix(""))
                            break
                elif path.is_dir():
                    add_candidate(path)
        except OSError:
            pass

        # Extracted NCBI archives often create one directory per volume:
        #   refseq_rna.00/refseq_rna.00.nhr
        #   refseq_rna.01/refseq_rna.01.nhr
        try:
            for child in container.iterdir():
                if not child.is_dir() or not volume_pattern.fullmatch(child.name):
                    continue
                volume_prefix = child / child.name
                if any(
                    Path(str(volume_prefix) + extension).is_file()
                    for extension in extensions
                ):
                    add_candidate(volume_prefix)
        except OSError:
            pass

        # Last-resort recursive discovery handles an additional nesting level
        # while remaining restricted to files belonging to this database name.
        try:
            for path in container.rglob("*"):
                if not path.is_file():
                    continue
                if not any(path.name.lower().endswith(ext.lower()) for ext in extensions):
                    continue
                stem = path.name
                for extension in extensions:
                    if stem.lower().endswith(extension.lower()):
                        stem = stem[: -len(extension)]
                        break
                candidate = path.with_name(stem)
                if volume_pattern.fullmatch(candidate.name):
                    add_candidate(candidate)
        except OSError:
            pass

        return dict(sorted(found.items(), key=lambda item: (int(item[0].rsplit('.', 1)[1]), item[0])))

    def _discover_files(
        self,
        name: str,
        database_type: str,
        *,
        prefix: Optional[Path] = None,
    ) -> List[Path]:
        prefix = prefix or self._find_existing_prefix(name, database_type)
        extensions = self._expected_extensions(database_type)
        files: List[Path] = []

        for extension in extensions:
            direct = Path(str(prefix) + extension)
            if direct.is_file():
                files.append(direct)

        volume_discovery_prefix = self._volume_container_prefix(prefix)
        for volume_prefix in self._discover_volume_prefixes(volume_discovery_prefix, database_type).values():
            for extension in extensions:
                candidate = Path(str(volume_prefix) + extension)
                if candidate.is_file():
                    files.append(candidate)

        unique = {path.resolve(): path for path in files}
        return sorted(unique.values(), key=lambda path: str(path).lower())

    def _resolve_expected_volume_prefix(
        self,
        alias_prefix: Path,
        volume_name: str,
        database_type: str,
    ) -> Optional[Path]:
        """Resolve one alias-declared BLAST volume deterministically.

        NCBI extracted archives used by Synapse store each volume in a
        directory such as ``refseq_rna.00/refseq_rna.00.nsq``.  Do not rely
        on recursive discovery for an alias-declared volume: resolve the
        exact volume name first, then accept the flat equivalent as a
        compatibility layout.
        """
        extensions = self._expected_extensions(database_type)
        container = alias_prefix.parent
        volume_dir = container / volume_name
        nested_prefix = volume_dir / volume_name

        if volume_dir.is_dir() and any(
            Path(str(nested_prefix) + extension).is_file()
            for extension in extensions
        ):
            return nested_prefix.resolve()

        flat_prefix = container / volume_name
        if any(
            Path(str(flat_prefix) + extension).is_file()
            for extension in extensions
        ):
            return flat_prefix.resolve()

        # Compatibility with an additional extraction directory level, but
        # only for the exact volume name declared by the alias.
        try:
            for candidate in container.rglob(volume_name):
                if not candidate.is_dir():
                    continue
                nested = candidate / volume_name
                if any(
                    Path(str(nested) + extension).is_file()
                    for extension in extensions
                ):
                    return nested.resolve()
        except OSError:
            pass

        return None


    def _missing_volume_files(
        self,
        volume_prefix: Path,
        database_type: str,
    ) -> List[str]:
        """Return missing mandatory sequence-volume files.

        NCBI multi-volume databases use the sequence core files for each
        volume. Version-5 metadata/index files are not necessarily present
        in every volume, so they are not required for volumes after the first.
        """
        if database_type == "protein":
            required = (".phr", ".pin", ".psq")
        else:
            required = (".nhr", ".nin", ".nsq")
        return [
            extension
            for extension in required
            if not Path(str(volume_prefix) + extension).is_file()
        ]

    def _volume_core_is_complete(
        self,
        volume_prefix: Path,
        database_type: str,
    ) -> bool:
        """Check the mandatory sequence-volume files for one volume."""
        return not self._missing_volume_files(volume_prefix, database_type)


    def scan_database(self, name: str, database_type: Optional[str] = None) -> LocalBlastDatabase:
        safe_name = self._safe_name(name)
        manifest = self._load_manifest()

        metadata = (
            manifest.get(safe_name, {})
            if isinstance(manifest.get(safe_name, {}), dict)
            else {}
        )
        catalog = LOCAL_DATABASE_CATALOG.get(safe_name, {})
        db_type = database_type or metadata.get("database_type") or catalog.get("type")

        if db_type not in {"nucleotide", "protein"}:
            direct_nucl = self._discover_files(safe_name, "nucleotide")
            direct_prot = self._discover_files(safe_name, "protein")
            if direct_prot and not direct_nucl:
                db_type = "protein"
            else:
                db_type = "nucleotide"

        prefix = self._find_existing_prefix(safe_name, db_type)
        files = self._discover_files(safe_name, db_type, prefix=prefix)

        alias_extension = ".pal" if db_type == "protein" else ".nal"
        alias = Path(str(prefix) + alias_extension)

        if db_type == "protein":
            direct_core = all(
                Path(str(prefix) + ext).is_file()
                for ext in (".phr", ".pin", ".psq")
            )
            v5_core = Path(str(prefix) + ".pdb").is_file()
        else:
            direct_core = all(
                Path(str(prefix) + ext).is_file()
                for ext in (".nhr", ".nin", ".nsq")
            )
            v5_core = Path(str(prefix) + ".ndb").is_file()

        volume_discovery_prefix = self._volume_container_prefix(prefix)
        discovered_volume_prefixes = self._discover_volume_prefixes(
            volume_discovery_prefix, db_type
        )
        expected_volume_names = _parse_blast_alias_volume_names(alias)

        # For an alias database, DBLIST is authoritative. Resolve each declared
        # volume by its exact name instead of trusting generic recursive
        # discovery. This is critical for extracted NCBI layouts such as:
        #   refseq_rna/refseq_rna.00/refseq_rna.00.nsq
        #   refseq_rna/refseq_rna.01/refseq_rna.01.nsq
        if expected_volume_names:
            volume_prefixes: Dict[str, Path] = {}
            for volume_name in expected_volume_names:
                resolved_volume = self._resolve_expected_volume_prefix(
                    volume_discovery_prefix, volume_name, db_type
                )

                # IMPORTANT: NCBI BLAST aliases may use an unpadded numeric
                # volume token (``.0``, ``.1`` ...) while the extracted archive
                # is stored with two digits (``.00``, ``.01`` ...). Treat both
                # representations as the same volume. Never construct a new
                # filesystem name by dropping zero padding. Resolve against
                # the actually discovered directory/file prefix instead.
                if resolved_volume is None:
                    expected_index = _extract_blast_volume_index(volume_name)
                    if expected_index is not None:
                        for observed_name, observed_prefix in discovered_volume_prefixes.items():
                            observed_index = _extract_blast_volume_index(observed_name)
                            if observed_index == expected_index:
                                resolved_volume = observed_prefix
                                break

                if resolved_volume is not None:
                    volume_prefixes[volume_name] = resolved_volume
        else:
            volume_prefixes = discovered_volume_prefixes

        observed_volume_names = list(volume_prefixes)

        missing_volumes = [
            volume
            for volume in expected_volume_names
            if volume not in volume_prefixes
        ]
        missing_volume_files: Dict[str, List[str]] = {}
        incomplete_volumes: List[str] = []
        for volume, volume_prefix in volume_prefixes.items():
            missing = self._missing_volume_files(volume_prefix, db_type)
            if missing:
                incomplete_volumes.append(volume)
                missing_volume_files[volume] = missing

        if expected_volume_names:
            all_expected_present = not missing_volumes
            all_expected_complete = all(
                volume in volume_prefixes
                and self._volume_core_is_complete(volume_prefixes[volume], db_type)
                for volume in expected_volume_names
            )
            structural = alias.is_file() and all_expected_present and all_expected_complete
        elif volume_prefixes:
            structural = alias.is_file() and not incomplete_volumes
        else:
            structural = alias.is_file() or direct_core or v5_core

        modified = None
        size = 0
        existing_files = [path for path in files if path.exists()]
        if existing_files:
            size = sum(path.stat().st_size for path in existing_files)
            newest = max(path.stat().st_mtime for path in existing_files)
            modified = datetime.fromtimestamp(newest, timezone.utc).replace(
                microsecond=0
            ).isoformat()

        # A catalog database such as refseq_rna may physically live in
        # blast_db/refseq_rna/refseq_rna. Keep the logical name stable while
        # retaining the physical prefix in ``path``. For explicitly nested
        # custom names (e.g. local_test/synapse_test), preserve the nested name.
        resolved_name = safe_name
        try:
            resolved_relative = prefix.relative_to(self.root)
            if "/" in safe_name:
                resolved_name = resolved_relative.as_posix()
            else:
                parent_name = prefix.parent.name
                is_numbered_volume_parent = bool(
                    re.fullmatch(rf"{re.escape(safe_name)}\.\d+", parent_name, re.IGNORECASE)
                )
                is_standard_catalog_prefix = (
                    prefix.parent.name == safe_name and prefix.name == safe_name
                )
                if not is_numbered_volume_parent and not is_standard_catalog_prefix:
                    if str(resolved_relative.parent) != ".":
                        resolved_name = resolved_relative.as_posix()
        except ValueError:
            pass

        resolved_metadata = manifest.get(resolved_name, {})
        if isinstance(resolved_metadata, dict) and resolved_metadata:
            metadata = resolved_metadata

        if expected_volume_names:
            volume_status = (
                f"{len(expected_volume_names) - len(missing_volumes)}/"
                f"{len(expected_volume_names)} volumes present"
            )
        elif volume_prefixes:
            volume_status = f"{len(volume_prefixes)} volume(s) detected"
        else:
            volume_status = "Single-volume database"

        if structural:
            verification_status = "STRUCTURALLY_VALID"
            verification_message = (
                f"BLAST database structure detected. {volume_status}. "
                "Use blastdbcmd for definitive verification."
            )
        elif existing_files:
            verification_status = "INCOMPLETE"
            details: List[str] = []
            if missing_volumes:
                details.append("missing: " + ", ".join(missing_volumes))
            if incomplete_volumes:
                details.append("incomplete: " + ", ".join(incomplete_volumes))
            verification_message = (
                "BLAST database files were found, but the database is incomplete. "
                + ("; ".join(details) if details else volume_status)
            )
        else:
            verification_status = "NOT_INSTALLED"
            verification_message = "Required BLAST database files were not detected."

        return LocalBlastDatabase(
            name=resolved_name,
            database_type=db_type,
            path=str(prefix),
            source=metadata.get("source") or catalog.get("source") or "Local",
            installed=bool(existing_files),
            structurally_valid=structural,
            verification_status=verification_status,
            verification_message=verification_message,
            file_count=len(existing_files),
            size_bytes=size,
            modified_at=modified,
            volumes=observed_volume_names,
            volume_paths=[
                str(path.parent if path.parent.name == path.name else path)
                for path in volume_prefixes.values()
            ],
            volume_count=len(volume_prefixes),
            expected_volume_count=(len(expected_volume_names) if expected_volume_names else None),
            missing_volumes=missing_volumes,
            incomplete_volumes=incomplete_volumes,
            missing_volume_files=missing_volume_files,
            manifest_metadata=metadata,
        )


    def list_databases(self) -> List[LocalBlastDatabase]:
        manifest = self._load_manifest()
        names = set(LOCAL_DATABASE_CATALOG)
        names.update(
            name for name in manifest
            if name != LOCAL_DB_MANIFEST_NAME and isinstance(name, str)
        )

        # Discover custom/nested databases, but NEVER promote a BLAST volume
        # directory (for example refseq_rna.00/) to a database of its own.
        # NCBI multi-volume databases are represented by one logical alias
        # (refseq_rna.nal) plus refseq_rna.00, refseq_rna.01, ... volumes.
        try:
            patterns = ("*.nal", "*.pal", "*.ndb", "*.pdb")
            volume_dir_pattern = re.compile(r"^.+\.\d+$", re.IGNORECASE)

            for pattern in patterns:
                for marker in self.root.rglob(pattern):
                    if not marker.is_file():
                        continue
                    try:
                        relative = marker.relative_to(self.root)
                    except ValueError:
                        continue

                    parts = relative.parts
                    # Synapse's own temporary views (synapse_blast_runtime_*,
                    # synapse_blast_verify_*) are not databases. One left
                    # behind by an interrupted search would otherwise be
                    # listed as a second copy of the database (2026-11-18).
                    if any(part.startswith("synapse_blast_") for part in parts[:-1]):
                        continue
                    # Ignore any marker physically stored inside a numbered
                    # volume directory. Such files belong to the parent
                    # database and must not create a second database entry.
                    if any(volume_dir_pattern.fullmatch(part) for part in parts[:-1]):
                        continue

                    stem = marker.stem
                    parent = marker.parent

                    # A volume file such as refseq_rna.00.ndb must also never
                    # become a standalone database.
                    if re.fullmatch(r".+\.\d+", stem):
                        continue

                    # Accept the supported container layout:
                    #   blast_db/name/name.nal
                    # as the logical database name ``name`` rather than the
                    # physical prefix ``name/name``. Flat databases and
                    # genuinely nested custom databases keep their relative path.
                    if parent.name == stem and parent != self.root:
                        try:
                            name = parent.relative_to(self.root).as_posix()
                        except ValueError:
                            name = stem
                    else:
                        name = relative.with_suffix("").as_posix()

                    if name:
                        names.add(name)
        except OSError:
            pass

        databases = [self.scan_database(name) for name in sorted(names, key=str.lower)]

        # Deduplicate aliases that resolve to the same physical prefix.
        # Prefer the shortest logical/catalog name (for example ``refseq_rna``)
        # over accidental nested aliases such as ``refseq_rna/refseq_rna.00/...``.
        unique: Dict[str, LocalBlastDatabase] = {}
        for database in databases:
            key = str(Path(database.path).resolve()).lower()
            current = unique.get(key)
            if current is None:
                unique[key] = database
                continue
            current_parts = len(Path(current.name).parts)
            candidate_parts = len(Path(database.name).parts)
            if candidate_parts < current_parts:
                unique[key] = database
            elif candidate_parts == current_parts and database.name.lower() < current.name.lower():
                unique[key] = database
        return sorted(unique.values(), key=lambda item: item.name.lower())


    def database_exists(self, name: str) -> bool:
        return self.scan_database(name).exists

    def _runtime_alias_parent(self, database: LocalBlastDatabase) -> Path:
        """Return a stable directory in which a runtime alias can resolve volumes.

        For nested NCBI layouts the physical volume prefixes are siblings below
        the database container, so the temporary alias must live in that same
        container (or a child of it). This keeps DBLIST entries relative and
        independent of Windows paths containing spaces.
        """
        database_path = Path(database.path).expanduser().resolve()
        volume_directory = database_path.parent
        if re.fullmatch(r".+\.\d+", volume_directory.name):
            return volume_directory.parent
        return volume_directory

    def verify_database(self, name: str, run_blastdbcmd: bool = True) -> LocalBlastDatabase:
        database = self.scan_database(name)
        if not database.structurally_valid:
            database.verification_status = "FAILED"
            details: List[str] = []
            if database.missing_volumes:
                details.append("missing volumes: " + ", ".join(database.missing_volumes))
            if database.missing_volume_files:
                file_details = "; ".join(
                    f"{volume} [{', '.join(files)}]"
                    for volume, files in database.missing_volume_files.items()
                )
                details.append("missing files: " + file_details)
            database.verification_message = (
                "Database structure is incomplete. "
                + ("; ".join(details) if details else "No complete database core was detected.")
            )
            return database

        if not run_blastdbcmd:
            database.verification_status = "STRUCTURALLY_VALID"
            database.verification_message = "File structure is valid."
            return database

        try:
            blastdbcmd = self.executable_manager.require("blastdbcmd")
        except BlastAnalysisError as exc:
            database.verification_status = "STRUCTURALLY_VALID_TOOL_MISSING"
            database.verification_message = (
                "File structure is valid, but blastdbcmd is unavailable. "
                f"{exc}"
            )
            return database

        runtime_directory: Optional[tempfile.TemporaryDirectory[str]] = None
        try:
            if _needs_runtime_blast_alias(database):
                runtime_parent = self._runtime_alias_parent(database)
                runtime_directory = tempfile.TemporaryDirectory(
                    prefix="synapse_blast_verify_",
                    dir=str(runtime_parent),
                )
                database_reference, database_working_directory = _create_runtime_blast_alias(
                    database, runtime_directory.name
                )
            else:
                database_reference, database_working_directory = (
                    _build_local_database_command_reference(database.path)
                )

            if _needs_runtime_blast_alias(database):
                env = _build_local_blast_environment(
                    self._runtime_alias_parent(database)
                )
            else:
                env = _build_local_blast_environment(database.path)
            command = [blastdbcmd, "-db", database_reference, "-info"]
            completed = subprocess.run(
                command,
                cwd=str(database_working_directory),
                env=env,
                capture_output=True,
                text=True,
                timeout=120,
                check=False,
            )
        except (OSError, subprocess.SubprocessError, BlastAnalysisError) as exc:
            database.verification_status = "VERIFICATION_ERROR"
            database.verification_message = str(exc)
            return database
        finally:
            if runtime_directory is not None:
                runtime_directory.cleanup()

        output = (completed.stdout or "") + "\r\n" + (completed.stderr or "")
        if completed.returncode == 0:
            parsed = _parse_blastdbcmd_info(output)
            database.title = parsed["title"]
            database.sequence_count = parsed["sequence_count"]
            database.total_length = parsed["total_length"]
            database.blastdb_version = parsed["blastdb_version"]
            database.database_date = parsed.get("database_date")
            database.volumes = parsed["volumes"]
            database.verification_status = "VERIFIED"
            database.verification_message = output.strip() or "blastdbcmd verification succeeded."
        else:
            database.verification_status = "FAILED"
            database.verification_message = output.strip() or (
                f"blastdbcmd exited with code {completed.returncode}."
            )
        return database

    def update_ncbi_database(
        self,
        name: str,
        *,
        timeout: Optional[float] = None,
    ) -> LocalBlastDatabase:
        if name not in LOCAL_DATABASE_CATALOG:
            raise BlastAnalysisError(
                f"'{name}' is not in the Synapse NCBI database catalog."
            )

        script = self.executable_manager.resolve_update_script()
        if not script:
            raise BlastAnalysisError(
                "update_blastdb.pl was not found. Install the NCBI BLAST+ package "
                "and make sure its bin directory is discoverable."
            )

        perl = self.executable_manager.resolve_perl()
        if not perl:
            raise BlastAnalysisError(
                "Perl was not found. NCBI's update_blastdb.pl requires Perl."
            )

        self.ensure_root()
        command = [
            perl,
            script,
            "--passive",
            "--decompress",
            name,
        ]

        try:
            completed = subprocess.run(
                command,
                cwd=str(self.root),
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise BlastAnalysisError(
                f"Database update timed out after {timeout} seconds."
            ) from exc
        except OSError as exc:
            raise BlastAnalysisError(f"Could not start database update: {exc}") from exc

        output = ((completed.stdout or "") + "\r\n" + (completed.stderr or "")).strip()
        if completed.returncode != 0:
            raise BlastAnalysisError(
                f"NCBI database update failed with exit code {completed.returncode}.\r\n"
                f"{output[-4000:]}"
            )

        manifest = self._load_manifest()
        manifest[name] = {
            "database_type": LOCAL_DATABASE_CATALOG[name]["type"],
            "source": LOCAL_DATABASE_CATALOG[name]["source"],
            "description": LOCAL_DATABASE_CATALOG[name]["description"],
            "last_update_command": "update_blastdb.pl --passive --decompress " + name,
            "updated_at_utc": _utc_now(),
        }
        self._save_manifest(manifest)

        database = self.scan_database(name)
        database.verification_message = output[-4000:] or "Database update completed."
        return database

    def create_custom_database(
        self,
        fasta_path: str | Path,
        name: str,
        database_type: str,
        *,
        title: Optional[str] = None,
        overwrite: bool = False,
    ) -> LocalBlastDatabase:
        database_type = database_type.lower().strip()
        if database_type not in {"nucleotide", "protein"}:
            raise BlastAnalysisError("database_type must be 'nucleotide' or 'protein'.")

        fasta = Path(fasta_path).expanduser().resolve()
        if not fasta.is_file():
            raise BlastAnalysisError(f"FASTA file not found: {fasta}")

        safe_name = self._safe_name(name)
        existing = self.scan_database(safe_name, database_type)
        if existing.installed and not overwrite:
            raise BlastAnalysisError(
                f"Database '{safe_name}' already exists. "
                "Choose another name or enable overwrite."
            )

        makeblastdb = self.executable_manager.require("makeblastdb")
        self.ensure_root()
        prefix = self._prefix(safe_name)
        dbtype_flag = "prot" if database_type == "protein" else "nucl"

        command = [
            makeblastdb,
            "-in",
            str(fasta),
            "-dbtype",
            dbtype_flag,
            "-parse_seqids",
            "-out",
            str(prefix),
            "-title",
            title or safe_name,
        ]

        try:
            completed = subprocess.run(
                command,
                cwd=str(self.root),
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise BlastAnalysisError(f"Could not start makeblastdb: {exc}") from exc

        output = ((completed.stdout or "") + "\r\n" + (completed.stderr or "")).strip()
        if completed.returncode != 0:
            raise BlastAnalysisError(
                f"makeblastdb failed with exit code {completed.returncode}.\r\n"
                f"{output[-4000:]}"
            )

        manifest = self._load_manifest()
        manifest[safe_name] = {
            "database_type": database_type,
            "source": "Custom FASTA",
            "description": title or f"Custom FASTA database: {fasta.name}",
            "fasta_source": str(fasta),
            "created_at_utc": _utc_now(),
        }
        self._save_manifest(manifest)

        database = self.scan_database(safe_name, database_type)
        database.verification_message = output[-4000:] or "Custom database created."
        return database

    def remove_database(self, name: str) -> List[str]:
        safe_name = self._safe_name(name)
        database = self.scan_database(safe_name)
        removed: List[str] = []

        if not database.installed:
            return removed

        prefix = Path(database.path).resolve()
        database_type = database.database_type
        files = self._discover_files(safe_name, database_type, prefix=prefix)
        volume_prefixes = self._discover_volume_prefixes(prefix, database_type)

        # Remove only artifacts belonging to this database. Archive files,
        # unrelated databases and the parent database container are preserved.
        targets: set[Path] = {path.resolve() for path in files if path.is_file()}
        targets.update(path.resolve() for path in volume_prefixes.values() if path.is_dir())

        # If the logical database prefix itself is a directory (unusual but
        # supported), remove only its database-owned files, never the parent.
        for target in sorted(targets, key=lambda path: (len(path.parts), str(path).lower()), reverse=True):
            try:
                if target.is_file():
                    target.unlink()
                    removed.append(str(target))
                elif target.is_dir():
                    target.rmdir()
                    removed.append(str(target))
            except OSError as exc:
                raise BlastAnalysisError(f"Could not remove '{target}': {exc}") from exc

        manifest = self._load_manifest()
        for manifest_name in {safe_name, database.name}:
            if manifest_name in manifest:
                del manifest[manifest_name]
        self._save_manifest(manifest)

        return removed


    def database_info(self, name: str) -> Dict[str, Any]:
        database = self.scan_database(name)
        return {
            "name": database.name,
            "database_type": database.database_type,
            "path": database.path,
            "source": database.source,
            "installed": database.installed,
            "structurally_valid": database.structurally_valid,
            "verification_status": database.verification_status,
            "verification_message": database.verification_message,
            "file_count": database.file_count,
            "size_bytes": database.size_bytes,
            "size_mb": round(database.size_mb, 2),
            "size_gb": round(database.size_gb, 3),
            "modified_at": database.modified_at,
            "title": database.title,
            "sequence_count": database.sequence_count,
            "total_length": database.total_length,
            "blastdb_version": database.blastdb_version,
            "volumes": database.volumes,
            "volume_paths": database.volume_paths,
            "volume_count": database.volume_count,
            "expected_volume_count": database.expected_volume_count,
            "missing_volumes": database.missing_volumes,
            "incomplete_volumes": database.incomplete_volumes,
            "missing_volume_files": database.missing_volume_files,
            "manifest_metadata": database.manifest_metadata,
        }
