"""Shared terminal presentation for Synapse (Stage F; audit D17).

Before this module every component drew its own frame: eleven different
separator widths were in use (42, 45, 58, 60, 65, 68, 70, 72, 78, 88, 90),
``[0]`` meant seven different things, and error messages appeared as
"ERROR:", "Error:" and "Error ". The modules looked like separate programs.

This module is the single source of that visual language. It holds no
scientific logic: it decides how a header, a menu, a prompt or a message is
drawn, never what it says. Each module keeps its own wording and its own
sections, which is what makes it recognizable as ORF analysis rather than
CRISPR design.

Conventions
-----------
Width
    One rule width, ``RULE_WIDTH``. 70 columns is the most common width in
    the existing code and leaves margin inside an 80-column terminal.
    ``MAJOR`` ("=") frames a screen, ``MINOR`` ("-") divides a section.

Navigation
    ``[0]`` always leaves the current screen. In a menu it is "Back"; in a
    dialog that would otherwise commit an action it is "Cancel"; only the
    main menu uses "Exit". No other meaning is given to ``[0]``, so a user
    can always press 0 to get out without losing work.

Messages
    One prefix per category, so a reader can tell at a glance whether
    something failed, needs attention or is merely informative:
    ``ERROR:``, ``WARNING:``, ``NOTE:``, ``OK:``.
"""

from typing import Iterable, Optional, Sequence, Tuple

# The string-building primitives live in core.textformat so that reports
# built by the analysis modules use the same layout without importing the
# terminal interface (P2-2). They are re-exported here unchanged.
from core.textformat import (  # noqa: F401
    FIELD_WIDTH,
    MAJOR,
    MINOR,
    RULE_WIDTH,
    columns,
    field_line,
    major_rule,
    minor_rule,
    rule,
)

BACK_KEY = "0"
BACK_LABEL = "Back"
CANCEL_LABEL = "Cancel"
EXIT_LABEL = "Exit"

ERROR_PREFIX = "ERROR:"
WARNING_PREFIX = "WARNING:"
NOTE_PREFIX = "NOTE:"
SUCCESS_PREFIX = "OK:"

CONTINUE_PROMPT = "Press Enter to continue..."
SELECTION_PROMPT = "Selection"
INVALID_SELECTION = "Invalid selection."


# ---------------------------------------------------------------------------
# Headers
# ---------------------------------------------------------------------------

def header(title: str, *, subtitle: Optional[str] = None, blank_before: bool = True) -> None:
    """Frame the top of a screen: rule, title, rule, optional subtitle."""
    if blank_before:
        print()
    print(major_rule())
    print(str(title).upper())
    print(major_rule())
    if subtitle:
        print(subtitle)


def section(title: str, *, blank_before: bool = True) -> None:
    """Start a block inside a screen: title over a minor rule."""
    if blank_before:
        print()
    print(str(title).upper())
    print(minor_rule())


def field(label: str, value: object, *, width: int = FIELD_WIDTH) -> None:
    """One aligned key/value line."""
    print(field_line(label, value, width=width))


# ---------------------------------------------------------------------------
# Menus and navigation
# ---------------------------------------------------------------------------

def option(key: object, label: str) -> None:
    """One selectable menu entry."""
    print(f"[{key}] {label}")


def options(entries: Iterable[Tuple[object, str]]) -> None:
    for key, label in entries:
        option(key, label)


def group(title: str) -> None:
    """A labeled group of entries inside a menu: a blank line, then the label."""
    print()
    print(str(title).upper())


def menu(groups: Sequence[Tuple[Optional[str], Sequence[Tuple[object, str]]]], *,
         zero: Optional[str] = BACK_LABEL) -> None:
    """The body of a menu screen, from the first entry to the closing rule.

    ``groups`` holds (label or None, entries) pairs; a module chooses its own
    groups and entries (Master Prompt 5.3), the layout is shared: each group
    after a blank line, then a blank line, the ``[0]`` entry and a minor rule.
    ``zero=None`` is for a screen that cannot be left (the sequence
    selector): it offers no ``[0]`` rather than giving 0 another meaning.
    The caller asks for the selection with ``ask_choice()``.
    """
    for title, entries in groups:
        if title:
            group(title)
        else:
            print()
        options(entries)
    print()
    if zero is not None:
        option(BACK_KEY, zero)
    print(minor_rule())


def back(label: str = BACK_LABEL) -> None:
    """The [0] entry that leaves the current screen."""
    option(BACK_KEY, label)


def cancel() -> None:
    """The [0] entry of a dialog that would otherwise commit an action."""
    option(BACK_KEY, CANCEL_LABEL)


def exit_option() -> None:
    """The [0] entry of the main menu only."""
    option(BACK_KEY, EXIT_LABEL)


def is_back(choice: str) -> bool:
    """True when the user asked to leave the current screen."""
    return str(choice).strip() == BACK_KEY


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def ask(label: str = SELECTION_PROMPT, *, default: Optional[str] = None) -> str:
    """Prompt for input. A default is shown in brackets and used when empty."""
    suffix = f" [{default}]" if default is not None else ""
    answer = input(f"{label}{suffix}: ").strip()
    return answer or (default or "")


def ask_choice(label: str = SELECTION_PROMPT) -> str:
    """Prompt for a menu selection."""
    return input(f"{label}: ").strip()


def pause(message: str = CONTINUE_PROMPT) -> None:
    input(f"\n{message}")


def confirm(question: str, *, default_yes: bool = False) -> bool:
    """Ask a yes/no question. Anything but an explicit yes means no."""
    hint = "[Y/n]" if default_yes else "[y/N]"
    answer = input(f"{question} {hint}: ").strip().lower()
    if not answer:
        return default_yes
    return answer in {"y", "yes", "j", "ja"}


# ---------------------------------------------------------------------------
# Messages
# ---------------------------------------------------------------------------

def error(message: str) -> None:
    print(f"{ERROR_PREFIX} {message}")


def warning(message: str) -> None:
    print(f"{WARNING_PREFIX} {message}")


def note(message: str) -> None:
    print(f"{NOTE_PREFIX} {message}")


def success(message: str) -> None:
    print(f"{SUCCESS_PREFIX} {message}")


def invalid_selection() -> None:
    print(INVALID_SELECTION)


def notes(messages: Sequence[str], *, prefix: str = NOTE_PREFIX) -> None:
    """Print a list of notes or warnings, one per line, or nothing if empty."""
    for message in messages:
        print(f"{prefix} {message}")
