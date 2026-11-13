"""NCBI BLAST polling survives a brief network failure and never loses the RID (D25).

Found on 2026-11-12 in a live search: NCBI accepted the query and returned
a request ID, then one status request failed with a DNS error
("getaddrinfo failed"). The first such error ended the search, and the
message did not name the RID, so a search that was still running at NCBI
could not be collected. No network is used here: the client's requests
are replaced.
"""

from __future__ import annotations

import logging
import unittest
from unittest import mock
from urllib.error import HTTPError, URLError

from analyses import blast as B


def _error(cause: BaseException) -> B.BlastAnalysisError:
    try:
        raise B.BlastAnalysisError(f"NCBI BLAST connection error: {cause}") from cause
    except B.BlastAnalysisError as error:
        return error


def _http(code: int, reason: str) -> HTTPError:
    error = HTTPError("https://blast.ncbi.nlm.nih.gov", code, reason, None, None)  # type: ignore[arg-type]
    error.close()
    return error


TRANSIENT = _error(URLError("[Errno 11001] getaddrinfo failed"))
SERVER = _error(_http(503, "busy"))
CLIENT = _error(_http(400, "bad request"))


class StatusPollingTests(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.CRITICAL)
        self.addCleanup(logging.disable, logging.NOTSET)
        self.client = B.NCBIBlastClient()
        patcher = mock.patch.object(B.time, "sleep", lambda seconds: None)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_a_brief_network_failure_is_retried(self):
        answers = [TRANSIENT, ("WAITING", ""), SERVER, ("READY", "")]

        def status(rid):
            answer = answers.pop(0)
            if isinstance(answer, Exception):
                raise answer
            return answer

        with mock.patch.object(self.client, "status", side_effect=status) as spy:
            self.client.await_ready("RID1", poll_interval=1, max_wait=600)
        self.assertEqual(spy.call_count, 4)

    def test_repeated_failures_end_the_wait(self):
        with mock.patch.object(self.client, "status", side_effect=TRANSIENT) as spy:
            with self.assertRaises(B.BlastAnalysisError):
                self.client.await_ready("RID1", poll_interval=1, max_wait=600)
        self.assertEqual(spy.call_count, B.STATUS_RETRIES + 1)

    def test_a_client_error_is_not_retried(self):
        with mock.patch.object(self.client, "status", side_effect=CLIENT) as spy:
            with self.assertRaises(B.BlastAnalysisError):
                self.client.await_ready("RID1", poll_interval=1, max_wait=600)
        self.assertEqual(spy.call_count, 1)

    def test_transient_errors_are_classified(self):
        self.assertTrue(B._is_transient_ncbi_error(TRANSIENT))
        self.assertTrue(B._is_transient_ncbi_error(SERVER))
        self.assertTrue(B._is_transient_ncbi_error(_error(TimeoutError())))
        self.assertTrue(B._is_transient_ncbi_error(_error(_http(429, "slow down"))))
        self.assertFalse(B._is_transient_ncbi_error(CLIENT))
        self.assertFalse(B._is_transient_ncbi_error(B.BlastAnalysisError("no cause")))


class SubmittedSearchTests(unittest.TestCase):
    def test_a_failed_search_names_its_request_id(self):
        logging.disable(logging.CRITICAL)
        self.addCleanup(logging.disable, logging.NOTSET)
        client = B.NCBIBlastClient()
        query = B.build_query("ACGT" * 15, name="q")
        with mock.patch.object(B.time, "sleep", lambda seconds: None), \
                mock.patch.object(client, "submit", return_value=("RIDXYZ123", 0, "")), \
                mock.patch.object(client, "status", side_effect=TRANSIENT):
            with self.assertRaises(B.BlastAnalysisError) as caught:
                client.run(query, B.BlastParameters(), poll_interval=1, max_wait=600)
        message = str(caught.exception)
        self.assertIn("RIDXYZ123", message)
        self.assertIn("Resume Search by Request ID", message)
        self.assertIn("getaddrinfo failed", message)


if __name__ == "__main__":
    unittest.main()
