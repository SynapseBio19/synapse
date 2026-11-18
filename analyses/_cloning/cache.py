"""Result cache for one strategy search (P2-6)."""

from __future__ import annotations

import copy
from contextlib import contextmanager
from typing import Any, Callable, TypeVar


# ============================================================
# STRATEGY SEARCH CACHE
# ============================================================

# The planner evaluates tens of thousands of strategies over the same two
# molecules, and every evaluation used to re-validate both sequences and
# re-scan them for each enzyme (P2-6). Both steps are pure functions of
# their arguments. Inside _strategy_search_cache() each distinct result is
# computed once and every caller receives its own deep copy, so no result
# shares objects with another. Outside the context nothing is cached.
_STRATEGY_SEARCH_CACHE = None


@contextmanager
def _strategy_search_cache():
    global _STRATEGY_SEARCH_CACHE
    if _STRATEGY_SEARCH_CACHE is not None:   # already inside a search
        yield
        return
    _STRATEGY_SEARCH_CACHE = {}
    try:
        yield
    finally:
        _STRATEGY_SEARCH_CACHE = None


_IMMUTABLE = (str, int, float, bool, type(None))
_T = TypeVar("_T")


def _copy_result(value: Any) -> Any:
    """An independent copy of a result tree (dicts, lists, plain values)."""
    if isinstance(value, _IMMUTABLE):
        return value
    if isinstance(value, dict):
        return {key: _copy_result(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_copy_result(item) for item in value]
    return copy.deepcopy(value)   # anything else: no shortcut


def _cached(key: tuple, compute: Callable[[], _T]) -> _T:
    cache = _STRATEGY_SEARCH_CACHE
    if cache is None:
        return compute()
    if key not in cache:
        cache[key] = compute()
    return _copy_result(cache[key])
