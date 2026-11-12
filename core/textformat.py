"""Plain-text layout primitives shared by reports and the terminal (P2-2).

Reports that the analysis modules build as strings (text exports, the
mutation and variant reports, the CRISPR report) and the terminal screens
use the same rule width and the same key/value alignment. These functions
only build strings; they print nothing and read no input, so the analysis
modules can use them without depending on the terminal interface.
``ui.style`` takes them from here, so the two can never drift apart.
"""

from typing import Sequence

RULE_WIDTH = 70

MAJOR = "="
MINOR = "-"

# Width of the label column in a key/value line, so aligned blocks line up
# across modules.
FIELD_WIDTH = 22


def rule(char: str = MAJOR, width: int = RULE_WIDTH) -> str:
    """One separator line. Prefer major_rule()/minor_rule()."""
    return char * width


def major_rule() -> str:
    return rule(MAJOR)


def minor_rule() -> str:
    return rule(MINOR)


def header_text(title: str) -> str:
    """The screen frame of ui.style.header as text: blank line, rule, title, rule."""
    return "\n".join(["", major_rule(), str(title).upper(), major_rule()])


def field_line(label: str, value: object, *, width: int = FIELD_WIDTH) -> str:
    """One aligned key/value line as text, for reports built as strings."""
    return f"{str(label) + ':':<{width}} {value}"


def columns(values: Sequence[object], widths: Sequence[int]) -> str:
    """One table row: each value left-aligned in its column, the last unpadded.

    Values are converted to text before padding, so a column stays aligned
    whatever the number of digits in it.
    """
    cells = [str(value) for value in values]
    padded = [f"{cell:<{width}}" for cell, width in zip(cells[:-1], widths)]
    return " ".join(padded + cells[len(padded):])
