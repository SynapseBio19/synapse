"""The HTTP layer: per-host throttling, bounded retry, JSON and URL helpers."""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Optional

from analyses._literature.models import DEFAULT_TIMEOUT, LiteratureSourceError, USER_AGENT


# Minimum interval between requests to the same host (seconds). NCBI
# E-utilities allow 3 requests/s without an API key and 10/s with one
# (SYNAPSE_NCBI_API_KEY); the other services ask for moderate, polite use.
MIN_REQUEST_INTERVAL_SECONDS = {
    "eutils.ncbi.nlm.nih.gov": 0.34,
    "www.ebi.ac.uk": 0.2,
    "api.openalex.org": 0.2,
    "api.crossref.org": 0.2,
}
RETRY_STATUS_CODES = {429, 500, 502, 503, 504}
MAX_REQUEST_ATTEMPTS = 3
_LAST_REQUEST_AT: dict[str, float] = {}


def _ncbi_api_key() -> Optional[str]:
    key = os.getenv("SYNAPSE_NCBI_API_KEY", "").strip()
    return key or None


def _throttle(host: str) -> None:
    interval = MIN_REQUEST_INTERVAL_SECONDS.get(host, 0.2)
    if host == "eutils.ncbi.nlm.nih.gov" and _ncbi_api_key():
        interval = 0.1
    last = _LAST_REQUEST_AT.get(host)
    if last is not None:
        wait = interval - (time.monotonic() - last)
        if wait > 0:
            time.sleep(wait)
    _LAST_REQUEST_AT[host] = time.monotonic()


def _retry_delay(exc: urllib.error.HTTPError, attempt: int) -> float:
    retry_after = exc.headers.get("Retry-After") if exc.headers else None
    if retry_after and retry_after.strip().isdigit():
        return min(30.0, float(retry_after))
    return min(30.0, 2.0 ** attempt)


def _request_text(url: str, timeout: int = DEFAULT_TIMEOUT, headers: Optional[dict[str, str]] = None) -> str:
    """GET a URL with per-host throttling and bounded retry on 429/5xx."""

    request_headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}

    if headers:

        request_headers.update(headers)

    host = urllib.parse.urlsplit(url).netloc

    for attempt in range(1, MAX_REQUEST_ATTEMPTS + 1):

        _throttle(host)

        request = urllib.request.Request(url, headers=request_headers)

        try:

            with urllib.request.urlopen(request, timeout=timeout) as response:

                charset = response.headers.get_content_charset() or "utf-8"

                return response.read().decode(charset, errors="replace")

        except urllib.error.HTTPError as exc:

            # The error carries the open response; its headers stay readable after close().
            exc.close()

            if exc.code in RETRY_STATUS_CODES and attempt < MAX_REQUEST_ATTEMPTS:

                time.sleep(_retry_delay(exc, attempt))

                continue

            raise LiteratureSourceError(

                f"HTTP {exc.code} from {host} after {attempt} attempt(s); "

                "the source may be rate-limiting requests or temporarily unavailable."

            ) from exc

        except (urllib.error.URLError, TimeoutError, OSError) as exc:

            raise LiteratureSourceError(f"Request failed for {host}: {exc}") from exc

    raise LiteratureSourceError(f"Request failed for {host}.")

def _request_json(url: str, timeout: int = DEFAULT_TIMEOUT, headers: Optional[dict[str, str]] = None) -> dict[str, Any]:

    try:

        payload = json.loads(_request_text(url, timeout=timeout, headers=headers))

    except json.JSONDecodeError as exc:

        raise LiteratureSourceError("The source returned invalid JSON.") from exc

    if not isinstance(payload, dict):

        raise LiteratureSourceError("The source returned an unexpected JSON structure.")

    return payload

def _url(base: str, params: dict[str, Any]) -> str:

    return base + "?" + urllib.parse.urlencode({k: v for k, v in params.items() if v is not None})


def _http_error_is_not_found(exc: Exception) -> bool:
    """Return True for a definite HTTP 404 response."""
    return "HTTP 404" in str(exc) or "status code 404" in str(exc).casefold()
