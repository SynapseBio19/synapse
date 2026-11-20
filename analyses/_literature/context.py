"""Literature context of a selected sequence or ORF."""

from __future__ import annotations

import re
from typing import Any, Optional


def build_sequence_literature_context(header: Optional[str] = None, orf: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    """Derive an optional search suggestion from the active sequence.

    Only the free-text sequence header (and the ORF label, if any) is used;
    no sequence data leaves the program. The suggestion is shown to the user
    and is never submitted without their confirmation.
    """
    context: dict[str, Any] = {}
    text = (header or "").strip()
    words = [w for w in re.split(r"[\s|,;_]+", text) if w]
    words = [w for w in words if not re.fullmatch(r"(?i)(topology|len|length)=.*|\d+", w)]
    if words:
        context["sequence_header"] = text
        context["suggested_query"] = " ".join(words[:6])
    if orf:
        context["orf"] = f"rank {orf.get('rank', '?')} frame {orf.get('frame', '?')}"
    return context
