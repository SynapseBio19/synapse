"""Client of the NCBI BLAST URL API: submission, polling with the NCBI rate
limits, result retrieval and resuming a saved RID."""

from __future__ import annotations

import logging
import re
import time
from typing import Any, Dict, Mapping, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from analyses._blast.config import (
    API_TIMEOUT_SECONDS,
    DEFAULT_EMAIL,
    DEFAULT_MAX_WAIT_SECONDS,
    DEFAULT_POLL_INTERVAL_SECONDS,
    DEFAULT_TOOL,
    NCBI_BLAST_URL,
    NCBI_MIN_POLL_INTERVAL_SECONDS,
    NCBI_MIN_REQUEST_INTERVAL_SECONDS,
)
from analyses._blast.models import (
    BlastAnalysisError,
    BlastParameters,
    BlastQuery,
    BlastRun,
    _utc_now,
)
from analyses._blast.query import check_program_query_type
from analyses._blast.parsing import describe_run_from_xml, parse_blast_xml


logger = logging.getLogger("synapse.blast")

# Consecutive failed status requests that are retried while a submitted
# search is polled (D25).
STATUS_RETRIES = 3


def _is_transient_ncbi_error(error: BaseException) -> bool:
    """True for a connection error, a timeout or an HTTP 429/5xx answer."""
    cause = error.__cause__
    if isinstance(cause, HTTPError):
        return cause.code == 429 or 500 <= cause.code < 600
    return isinstance(cause, (URLError, TimeoutError, OSError))


# ============================================================================
# NCBI BLAST API
# ============================================================================

class NCBIBlastClient:
    """
    Minimal dependency-free client for the NCBI BLAST URL API.

    The client performs:
        PUT -> receive RID
        GET -> poll status
        GET -> retrieve XML

    It intentionally does not silently fall back to fake/local hits.
    """

    def __init__(
        self,
        *,
        tool: str = DEFAULT_TOOL,
        email: str = DEFAULT_EMAIL,
        timeout: float = API_TIMEOUT_SECONDS,
        user_agent: Optional[str] = None,
        min_request_interval: float = NCBI_MIN_REQUEST_INTERVAL_SECONDS,
        min_poll_interval: float = NCBI_MIN_POLL_INTERVAL_SECONDS,
    ) -> None:
        raw_tool = tool.strip() or DEFAULT_TOOL
        # NCBI recommends a tool identifier without internal spaces.
        self.tool = re.sub(r"\s+", "_", raw_tool)
        self.email = email.strip()
        self.timeout = float(timeout)
        self.min_request_interval = max(0.0, float(min_request_interval))
        self.min_poll_interval = max(1.0, float(min_poll_interval))
        self.user_agent = user_agent or f"{self.tool}/1.0"
        self._last_request_at: Optional[float] = None

    def _throttle_request(self) -> None:
        """Respect NCBI's public BLAST request-frequency guidance."""
        if self._last_request_at is None:
            return
        elapsed = time.monotonic() - self._last_request_at
        remaining = self.min_request_interval - elapsed
        if remaining > 0:
            time.sleep(remaining)

    def _request(
        self,
        params: Mapping[str, Any],
        *,
        method: str = "GET",
    ) -> str:
        encoded = urlencode(
            [(key, value) for key, value in params.items() if value is not None]
        )

        if method.upper() == "POST":
            request = Request(
                NCBI_BLAST_URL,
                data=encoded.encode("utf-8"),
                method="POST",
            )
        else:
            separator = "&" if "?" in NCBI_BLAST_URL else "?"
            request = Request(
                f"{NCBI_BLAST_URL}{separator}{encoded}",
                method="GET",
            )

        request.add_header("User-Agent", self.user_agent)
        request.add_header("Accept", "text/plain, application/xml, */*")

        self._throttle_request()
        try:
            with urlopen(request, timeout=self.timeout) as response:
                self._last_request_at = time.monotonic()
                return response.read().decode("utf-8", errors="replace")
        except HTTPError as exc:
            try:
                body = exc.read().decode("utf-8", errors="replace")
            finally:
                exc.close()
            raise BlastAnalysisError(
                f"NCBI BLAST HTTP error {exc.code}: {body[:500]}"
            ) from exc
        except URLError as exc:
            raise BlastAnalysisError(f"NCBI BLAST connection error: {exc}") from exc
        except TimeoutError as exc:
            raise BlastAnalysisError("NCBI BLAST request timed out.") from exc

    def submit(
        self,
        query: BlastQuery,
        parameters: BlastParameters,
    ) -> Tuple[str, Optional[float], str]:
        params = parameters.normalized()
        payload: Dict[str, Any] = {
            "CMD": "Put",
            "PROGRAM": params.program,
            "DATABASE": params.database,
            "QUERY": f">{query.name}\r\n{query.sequence}",
            "EXPECT": params.evalue,
            "HITLIST_SIZE": params.max_target_seqs,
            "TOOL": params.tool,
        }

        if params.email:
            payload["EMAIL"] = params.email
        if params.word_size is not None:
            payload["WORD_SIZE"] = params.word_size
        if params.gap_open is not None:
            payload["GAPOPEN"] = params.gap_open
        if params.gap_extend is not None:
            payload["GAPEXT"] = params.gap_extend
        if params.matrix:
            payload["MATRIX_NAME"] = params.matrix
        if params.filter_query:
            payload["FILTER"] = params.filter_query

        response = self._request(payload, method="POST")

        rid_match = re.search(r"RID\s*=\s*([A-Z0-9_-]+)", response, re.IGNORECASE)
        if not rid_match:
            raise BlastAnalysisError(
                "NCBI did not return a request ID (RID). "
                f"Response: {response[:1000]}"
            )

        rtoe_match = re.search(r"RTOE\s*=\s*(\d+)", response, re.IGNORECASE)
        rtoe = float(rtoe_match.group(1)) if rtoe_match else None
        rid = rid_match.group(1)
        return rid, rtoe, response

    def status(self, rid: str) -> Tuple[str, str]:
        response = self._request(
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo",
            }
        )

        status_match = re.search(
            r"Status\s*=\s*([A-Za-z]+)", response, re.IGNORECASE
        )
        status = status_match.group(1).upper() if status_match else "UNKNOWN"
        return status, response

    def fetch_xml(self, rid: str) -> str:
        response = self._request(
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "XML",
            }
        )
        if "<BlastOutput" not in response and "BlastOutput2" not in response:
            raise BlastAnalysisError(
                "NCBI returned no recognizable BLAST XML result. "
                f"Response: {response[:1000]}"
            )
        return response

    def run(
        self,
        query: BlastQuery,
        parameters: BlastParameters,
        *,
        poll_interval: float = DEFAULT_POLL_INTERVAL_SECONDS,
        max_wait: float = DEFAULT_MAX_WAIT_SECONDS,
    ) -> BlastRun:
        params = parameters.normalized()
        check_program_query_type(params.program, query)

        if params.email == "":
            # NCBI can accept requests without email, but Synapse exposes this
            # explicitly instead of inventing one.
            email_warning = (
                "No contact email was configured for the NCBI request. "
                "Set an email address in Synapse configuration when appropriate."
            )
        else:
            email_warning = None

        submitted_at = _utc_now()
        rid, rtoe, _submit_response = self.submit(query, params)

        # NCBI recommends no more than one status request per RID per minute.
        # RTOE is only an estimate, so we honor it while enforcing the minimum
        # polling interval as the hard lower bound.
        wait_time = max(0.0, rtoe or 0.0, self.min_poll_interval)
        if wait_time:
            time.sleep(min(wait_time, max_wait))

        try:
            run = self.resume(
                rid,
                query=query,
                parameters=params,
                poll_interval=poll_interval,
                max_wait=max_wait,
                submitted_at=submitted_at,
            )
        except BlastAnalysisError as exc:
            # The search exists at NCBI once it has an RID; name it, so that it
            # can be collected later instead of being submitted again (D25).
            raise BlastAnalysisError(
                f"{exc} The search was submitted to NCBI as request ID {rid}; "
                "collect it later with Resume Search by Request ID (BLAST menu [2]) "
                "or core.api.blast_resume while NCBI keeps the RID."
            ) from exc

        if email_warning:
            run.metadata.warnings.append(email_warning)

        return run

    def await_ready(
        self,
        rid: str,
        *,
        poll_interval: float = DEFAULT_POLL_INTERVAL_SECONDS,
        max_wait: float = DEFAULT_MAX_WAIT_SECONDS,
    ) -> None:
        """Poll a request ID until NCBI reports READY, or raise.

        A connection error, a timeout or an HTTP 429/5xx answer to a status
        request is retried up to ``STATUS_RETRIES`` times in a row, within
        ``max_wait``; the search keeps running at NCBI meanwhile. Before
        2026-11-12 the first such error ended the search (D25).
        """
        start = time.monotonic()
        failures = 0

        while True:
            try:
                status, _status_response = self.status(rid)
                failures = 0
            except BlastAnalysisError as exc:
                if not _is_transient_ncbi_error(exc) or failures >= STATUS_RETRIES \
                        or time.monotonic() - start >= max_wait:
                    raise
                failures += 1
                logger.warning("NCBI status request for %s failed (%s); retry %d of %d.",
                               rid, exc, failures, STATUS_RETRIES)
                time.sleep(max(self.min_poll_interval, poll_interval))
                continue

            if status == "READY":
                return
            if status in {"FAILED", "UNKNOWN"}:
                raise BlastAnalysisError(
                    f"NCBI BLAST request {rid} ended with status {status}. "
                    "NCBI keeps a request ID only for a limited time, so an "
                    "older RID may simply have expired."
                )
            if time.monotonic() - start >= max_wait:
                raise BlastAnalysisError(
                    f"BLAST request {rid} exceeded the maximum wait time "
                    f"of {max_wait:.0f} seconds."
                )

            time.sleep(max(self.min_poll_interval, poll_interval))

    def resume(
        self,
        rid: str,
        *,
        query: Optional[BlastQuery] = None,
        parameters: Optional[BlastParameters] = None,
        poll_interval: float = DEFAULT_POLL_INTERVAL_SECONDS,
        max_wait: float = DEFAULT_MAX_WAIT_SECONDS,
        submitted_at: Optional[str] = None,
    ) -> BlastRun:
        """Collect a search that was submitted earlier, by its request ID (P1-7).

        Every Synapse run records its RID, so a search can be picked up again
        without resubmitting the query - useful when a long search outlives
        the session that started it. NCBI retains an RID only for a limited
        period; an expired one reports UNKNOWN.

        ``query`` and ``parameters`` are used when the caller still holds
        them, for example from a stored run. Otherwise they are rebuilt from
        the returned XML, which states the query definition, the query
        length, the program and the database, but not the query sequence.
        """
        rid = str(rid).strip()
        if not rid:
            raise BlastAnalysisError("A request ID is required to resume a search.")

        self.await_ready(rid, poll_interval=poll_interval, max_wait=max_wait)
        xml_text = self.fetch_xml(rid)
        completed_at = _utc_now()

        if query is None or parameters is None:
            recovered_query, recovered_parameters = describe_run_from_xml(xml_text)
            query = query or recovered_query
            parameters = parameters or recovered_parameters

        run = parse_blast_xml(
            xml_text,
            query=query,
            rid=rid,
            parameters=parameters,
            submitted_at=submitted_at or completed_at,
            completed_at=completed_at,
        )
        if submitted_at is None:
            run.metadata.warnings.append(
                "This run was resumed from its request ID; the original "
                "submission time is not known to this session."
            )
        return run
