"""Size, version and disk-space checks for NCBI's preformatted BLAST databases.

NCBI publishes a small JSON file per database next to the archives
(``https://ftp.ncbi.nlm.nih.gov/blast/db/<name>-<nucl|prot>-metadata.json``,
BLAST database metadata version 1.1). It gives the number of sequences, the
date of the release, the download size (``bytes-total-compressed``) and the
size of the decompressed database (``bytes-total``). Synapse reads it before
a download so that the user sees what will be written and whether the disk
can hold it, and to tell whether NCBI has a newer release than the local
copy. Added on 2026-11-19 (Phase 9).
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any, Callable, Mapping, Optional
from urllib.error import URLError
from urllib.request import Request, urlopen

from analyses._blast.config import API_TIMEOUT_SECONDS, DEFAULT_TOOL, LOCAL_DATABASE_CATALOG
from analyses._blast.models import BlastAnalysisError, LocalBlastDatabase

NCBI_DB_METADATA_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/{name}-{kind}-metadata.json"
_KIND = {"nucleotide": "nucl", "protein": "prot"}


@dataclass(frozen=True)
class NcbiDatabaseMetadata:
    """What NCBI states about one release of a preformatted database."""

    name: str
    database_type: str
    description: str
    sequence_count: int
    letter_count: int
    last_updated: date
    download_bytes: int
    disk_bytes: int
    volume_count: int

    def required_free_bytes(self, installed_bytes: int = 0) -> int:
        """Free space needed for a download or update.

        The decompressed database (less what an installed copy already
        occupies, since an update replaces it) plus one volume archive:
        update_blastdb.pl fetches and unpacks one volume at a time and
        removes the archive afterwards. An estimate from NCBI's totals; the
        largest volume can be bigger than the average.
        """
        per_volume = math.ceil(self.download_bytes / max(1, self.volume_count))
        return max(0, self.disk_bytes - max(0, installed_bytes)) + per_volume


@dataclass(frozen=True)
class DownloadSpaceCheck:
    required_bytes: int
    free_bytes: int

    @property
    def sufficient(self) -> bool:
        return self.free_bytes >= self.required_bytes


def metadata_url(name: str) -> str:
    if name not in LOCAL_DATABASE_CATALOG:
        raise BlastAnalysisError(f"'{name}' is not in the Synapse NCBI database catalog.")
    kind = _KIND[LOCAL_DATABASE_CATALOG[name]["type"]]
    return NCBI_DB_METADATA_URL.format(name=name, kind=kind)


def parse_database_metadata(payload: Mapping[str, Any], name: str) -> NcbiDatabaseMetadata:
    """Read NCBI's metadata JSON; raises BlastAnalysisError if a field is missing."""
    try:
        updated = datetime.strptime(str(payload["last-updated"])[:10], "%Y-%m-%d").date()
        metadata = NcbiDatabaseMetadata(
            name=str(payload.get("dbname") or name),
            database_type=str(payload["dbtype"]).lower(),
            description=str(payload.get("description") or ""),
            sequence_count=int(payload["number-of-sequences"]),
            letter_count=int(payload["number-of-letters"]),
            last_updated=updated,
            download_bytes=int(payload["bytes-total-compressed"]),
            disk_bytes=int(payload["bytes-total"]),
            volume_count=int(payload["number-of-volumes"]),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise BlastAnalysisError(f"NCBI database metadata for '{name}' is incomplete: {exc}") from exc
    if metadata.name != name:
        raise BlastAnalysisError(f"NCBI returned metadata for '{metadata.name}', not '{name}'.")
    return metadata


def fetch_database_metadata(
    name: str,
    *,
    timeout: float = API_TIMEOUT_SECONDS,
    opener: Callable[..., Any] = urlopen,
) -> NcbiDatabaseMetadata:
    """Download and read NCBI's metadata file for ``name`` (a few kilobytes)."""
    request = Request(metadata_url(name), headers={"User-Agent": f"{DEFAULT_TOOL} (BLAST database metadata)"})
    try:
        with opener(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except (URLError, OSError, ValueError) as exc:
        raise BlastAnalysisError(f"Could not read NCBI's metadata for '{name}': {exc}") from exc
    return parse_database_metadata(payload, name)


def free_disk_bytes(path: str | Path) -> int:
    """Free space on the drive that holds ``path`` (or its nearest existing parent)."""
    candidate = Path(path).expanduser().resolve()
    while not candidate.exists() and candidate.parent != candidate:
        candidate = candidate.parent
    return int(shutil.disk_usage(candidate).free)


def check_download_space(metadata: NcbiDatabaseMetadata, root: str | Path,
                         installed_bytes: int = 0) -> DownloadSpaceCheck:
    return DownloadSpaceCheck(required_bytes=metadata.required_free_bytes(installed_bytes),
                              free_bytes=free_disk_bytes(root))


def local_release_date(database: LocalBlastDatabase) -> Optional[date]:
    """The release date blastdbcmd reported for a local database, if any.

    blastdbcmd -info prints it like ``Sep 30, 2026  8:36 PM``.
    """
    text = " ".join(str(database.database_date or "").split())
    match = re.match(r"([A-Z][a-z]{2}) (\d{1,2}), (\d{4})", text)
    if not match:
        return None
    try:
        return datetime.strptime(" ".join(match.groups()), "%b %d %Y").date()
    except ValueError:
        return None


def newer_release_available(database: LocalBlastDatabase, metadata: NcbiDatabaseMetadata) -> Optional[bool]:
    """True if NCBI's release is newer than the local copy; None if unknown."""
    local = local_release_date(database)
    if local is None:
        return None
    return metadata.last_updated > local


def format_bytes(value: int) -> str:
    return f"{value / 1e9:,.1f} GB" if value >= 1e8 else f"{value / 1e6:,.1f} MB"
