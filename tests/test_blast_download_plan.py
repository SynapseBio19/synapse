"""Size, release and disk-space checks before an NCBI database download (2026-11-19).

The fixture ``tests/fixtures/ncbi_blastdb/swissprot-prot-metadata.json`` is
NCBI's metadata file for swissprot as published on 2026-11-19 (see
SOURCE.md). Nothing is downloaded by these tests; the network is replaced
by a fake opener and the disk by a fake free-space value.
"""

from __future__ import annotations

from datetime import date
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
from urllib.error import URLError

from analyses import blast as B
from analyses._blast import ncbi_metadata
from ui import blast_menu

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "ncbi_blastdb" / "swissprot-prot-metadata.json"
PAYLOAD = json.loads(FIXTURE.read_text(encoding="utf-8"))


class _Opener:
    """Stands in for urlopen and remembers the requested URL."""

    def __init__(self, payload=None, error=None):
        self.payload, self.error, self.url = payload, error, ""

    def __call__(self, request, timeout):
        if self.error:
            raise self.error
        self.url = request.full_url
        return io.BytesIO(json.dumps(self.payload).encode("utf-8"))


def _opener(payload=None, error=None):
    return _Opener(payload, error)


class MetadataTests(unittest.TestCase):
    def test_ncbi_file_is_read(self):
        metadata = B.parse_database_metadata(PAYLOAD, "swissprot")
        self.assertEqual(metadata.database_type, "protein")
        self.assertEqual(metadata.sequence_count, 487942)
        self.assertEqual(metadata.last_updated, date(2026, 11, 8))
        self.assertEqual((metadata.download_bytes, metadata.disk_bytes, metadata.volume_count),
                         (225147306, 366074859, 1))
        self.assertEqual(metadata.required_free_bytes(), 366074859 + 225147306)
        self.assertEqual(metadata.required_free_bytes(installed_bytes=366074859), 225147306)

    def test_fetch_uses_the_catalog_type(self):
        opener = _opener(PAYLOAD)
        metadata = B.fetch_database_metadata("swissprot", opener=opener)
        self.assertEqual(opener.url, "https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot-prot-metadata.json")
        self.assertEqual(metadata.name, "swissprot")
        self.assertEqual(B.metadata_url("refseq_rna"),
                         "https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_rna-nucl-metadata.json")

    def test_failures_are_reported(self):
        with self.assertRaisesRegex(B.BlastAnalysisError, "Could not read NCBI's metadata for 'swissprot'"):
            B.fetch_database_metadata("swissprot", opener=_opener(error=URLError("offline")))
        incomplete = {k: v for k, v in PAYLOAD.items() if k != "bytes-total"}
        with self.assertRaisesRegex(B.BlastAnalysisError, "incomplete"):
            B.parse_database_metadata(incomplete, "swissprot")
        with self.assertRaisesRegex(B.BlastAnalysisError, "not 'nr'"):
            B.parse_database_metadata(PAYLOAD, "nr")
        with self.assertRaisesRegex(B.BlastAnalysisError, "not in the Synapse NCBI database catalog"):
            B.metadata_url("my_private_db")

    def test_release_comparison(self):
        metadata = B.parse_database_metadata(PAYLOAD, "swissprot")
        older = B.LocalBlastDatabase(name="swissprot", database_type="protein", path="x",
                                     database_date="Sep 30, 2026  8:36 PM")
        current = B.LocalBlastDatabase(name="swissprot", database_type="protein", path="x",
                                       database_date="Nov 8, 2026  1:00 AM")
        unknown = B.LocalBlastDatabase(name="swissprot", database_type="protein", path="x")
        self.assertEqual(B.local_release_date(older), date(2026, 9, 30))
        self.assertIs(B.newer_release_available(older, metadata), True)
        self.assertIs(B.newer_release_available(current, metadata), False)
        self.assertIsNone(B.newer_release_available(unknown, metadata))

    def test_free_space_of_a_folder_that_does_not_exist_yet(self):
        missing = Path(tempfile.gettempdir()) / "synapse_no_such_dir" / "deeper"
        self.assertGreater(B.free_disk_bytes(missing), 0)


class DownloadGuardTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="synapse_dl_"))
        self.addCleanup(lambda: __import__("shutil").rmtree(self.root, ignore_errors=True))
        self.manager = B.LocalBlastDatabaseManager(self.root)

    def test_update_refuses_to_start_without_space(self):
        with mock.patch.object(self.manager.executable_manager, "resolve_update_script", return_value="u.pl"), \
                mock.patch.object(self.manager.executable_manager, "resolve_perl", return_value="perl"), \
                mock.patch.object(ncbi_metadata, "free_disk_bytes", return_value=10 ** 9), \
                mock.patch("analyses._blast.databases.subprocess.run") as run:
            with self.assertRaisesRegex(B.BlastAnalysisError, "Not enough free disk space for refseq_rna"):
                self.manager.update_ncbi_database("refseq_rna", required_free_bytes=81 * 10 ** 9)
        run.assert_not_called()

    def _plan(self, free, payload=PAYLOAD, error=None):
        def fetch(name):
            if error:
                raise B.BlastAnalysisError(error)
            return B.parse_database_metadata(payload, name)
        with mock.patch.object(blast_menu, "fetch_database_metadata", side_effect=fetch), \
                mock.patch.object(ncbi_metadata, "free_disk_bytes", return_value=free), \
                mock.patch("builtins.print"):
            return blast_menu._show_ncbi_download_plan(self.manager, "swissprot")

    def test_the_dialog_offers_the_download_only_with_enough_space(self):
        self.assertEqual(self._plan(10 ** 12), 366074859 + 225147306)
        self.assertIs(self._plan(10 ** 6), False)

    def test_without_metadata_the_dialog_still_offers_the_download(self):
        self.assertIsNone(self._plan(10 ** 12, error="offline"))


if __name__ == "__main__":
    unittest.main()
