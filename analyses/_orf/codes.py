"""Start and stop codons, search policies, genetic-code tables and codon classification."""

from __future__ import annotations

from core.seqtools import codon_dict
from functools import lru_cache
from itertools import product


# ============================================================
# ORF ANALYSIS
# ============================================================


# ============================================================
# CONSTANTS
# ============================================================

STANDARD_START_CODONS = {"ATG"}

ALTERNATIVE_START_CODONS = {
    "GTG",
    "TTG"
}

ALL_START_CODONS = (
    STANDARD_START_CODONS
    | ALTERNATIVE_START_CODONS
)

STANDARD_STOP_CODONS = {
    "TAA",
    "TAG",
    "TGA"
}

# ------------------------------------------------------------
# Default search policy (audit P0-9, 2026-10-29)
# ------------------------------------------------------------
# Defaults follow NCBI ORFfinder conventions: minimal ORF length 75 nt
# (stop codon included), ATG-only start codons, and one ORF per stop codon
# (the longest, i.e. outermost in-frame start). These are screening
# conventions, not biological thresholds: short ORFs, non-ATG starts and
# nested ORFs exist and can be enabled explicitly.
#
# Previous defaults (still selectable): min length 10 nt, starts
# ATG/GTG/TTG, every in-frame start reported, heuristic score ranking.
DEFAULT_MIN_ORF_LENGTH = 75

START_CODON_POLICIES = {
    "atg_only": frozenset(STANDARD_START_CODONS),
    "atg_and_alternative": frozenset(ALL_START_CODONS),
}

DEFAULT_START_CODON_POLICY = "atg_only"
DEFAULT_REPORT_NESTED = False
DEFAULT_RANKING = "length"
ORF_RANKING_MODES = ("length", "heuristic")

LEGACY_ORF_SETTINGS = {
    "min_length": 10,
    "start_codon_policy": "atg_and_alternative",
    "nested": True,
    "ranking": "heuristic",
}


DNA_IUPAC = {
    "A": ("A",),
    "T": ("T",),
    "C": ("C",),
    "G": ("G",),

    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),

    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "T", "C", "G")
}


# ============================================================
# GENETIC CODE
# ============================================================

# Standard code (NCBI table 1) from core.seqtools; stops are "STOP" here.
GENETIC_CODE = codon_dict(1, stop_symbol="STOP")


# ============================================================
# IUPAC CODON HANDLING
# ============================================================

def expand_iupac_codon(codon):

    # A fresh list per call: callers may change it; the expansion itself
    # is cached (a genome-wide ORF search classifies millions of codons
    # drawn from a few hundred distinct strings).
    return list(_expanded_codon(codon))


@lru_cache(maxsize=4096)
def _expanded_codon(codon):

    codon = codon.upper()

    if len(codon) != 3:
        return ()

    possibilities = []

    for base in codon:

        if base not in DNA_IUPAC:
            return ()

        possibilities.append(
            DNA_IUPAC[base]
        )

    return tuple(
        "".join(combination)
        for combination in product(
            *possibilities
        )
    )


def _stop_codons_of(codon_table):
    if codon_table is None:
        return STANDARD_STOP_CODONS
    return {codon for codon, aa in codon_table.items() if aa == "STOP"}


def codon_table_for(genetic_code=1):
    """Return the local codon->amino-acid dict ('STOP' for stops) for an NCBI table.

    Table 1 returns the module's own GENETIC_CODE (cross-checked against
    core.seqtools in the test suite); other tables come from core.seqtools.
    """
    if genetic_code in (None, 1):
        return GENETIC_CODE
    from core.seqtools import get_genetic_code
    table = get_genetic_code(int(genetic_code))
    # Codons that can terminate count as stops, also where NCBI lists an
    # amino acid for them (context-dependent codons of codes 27, 28, 31).
    return {codon: ("STOP" if aa == "*" or codon in table.stop_codons else aa)
            for codon, aa in table.codons.items()}


def start_codons_for(genetic_code=1, start_codon_policy=None):
    """Start codons for a genetic code and a start-codon policy.

    Table 1 keeps the historical Synapse policies (ATG, or ATG+GTG+TTG).
    For any table, "atg_only" means ATG. For tables other than 1,
    "atg_and_alternative" uses the start codons NCBI lists for that table;
    "table_starts" does so for every table (table 1: TTG, CTG, ATG).
    """
    policy = start_codon_policy or DEFAULT_START_CODON_POLICY
    if policy not in START_CODON_POLICIES and policy != "table_starts":
        raise ValueError(
            f"Unknown start codon policy '{policy}'. "
            f"Use one of: {', '.join(list(START_CODON_POLICIES) + ['table_starts'])}."
        )
    if policy == "atg_only":
        return set(STANDARD_START_CODONS)
    if policy == "table_starts" or genetic_code not in (None, 1):
        from core.seqtools import get_genetic_code
        return set(get_genetic_code(int(genetic_code or 1)).start_codons)
    return set(START_CODON_POLICIES[policy])


def initiation_codons(genetic_code=1):
    """Codons that can initiate an ORF in a genetic code (audit D1).

    The union of the start codons of every start-codon policy for the
    table: the NCBI start codons of that table (core.seqtools, checked
    against Biopython) plus, for table 1, Synapse's documented
    "atg_and_alternative" policy (ATG, GTG, TTG). As the first codon of an
    ORF such a codon is the initiation codon and is decoded as methionine;
    anywhere else it is read with the elongation code.
    """
    policies = list(START_CODON_POLICIES) + ["table_starts"]
    return frozenset().union(
        *(start_codons_for(genetic_code, policy) for policy in policies)
    )


def is_initiation_codon(codon, genetic_code=1):
    """True when every IUPAC reading of the codon is an initiation codon.

    This is the rule the ORF finder applies to a "definite" start codon,
    so the codon that starts a Synapse ORF is always decoded as initiator.
    Accepts DNA or RNA letters.
    """
    codon = str(codon).upper().replace("U", "T")
    return analyze_codon_identity(codon, initiation_codons(genetic_code)) == "definite"


def _genetic_code_name(genetic_code=1):
    if genetic_code in (None, 1):
        return "NCBI translation table 1 (standard)"
    from core.seqtools import get_genetic_code
    table = get_genetic_code(int(genetic_code))
    return f"NCBI translation table {table.table_id} ({table.name})"


def analyze_codon_identity(
    codon,
    target_codons
):

    expanded = _expanded_codon(
        codon
    )

    if not expanded:
        return "invalid"

    matches = [
        expanded_codon in target_codons
        for expanded_codon in expanded
    ]

    if all(matches):
        return "definite"

    if any(matches):
        return "possible"

    return "no"


# ============================================================
# START CODON CLASSIFICATION
# ============================================================

def _is_alternative_start(orf):
    """A non-ATG start that initiates in the ORF's genetic code (audit S19)."""
    codon = orf["start_codon"]
    return codon != "ATG" and is_initiation_codon(codon, orf.get("genetic_code") or 1)


def classify_start_codon(codon, genetic_code=1):
    """Canonical (ATG), Alternative, Possible or None for a genetic code.

    The start set is ``initiation_codons(genetic_code)`` (audit S19; it was
    the fixed ATG/GTG/TTG set, so CTG in table 1 or ATA/ATT in table 2
    were "None").
    """

    status = analyze_codon_identity(
        codon,
        initiation_codons(genetic_code)
    )

    if status == "definite":

        if codon == "ATG":
            return "Canonical"

        if codon in initiation_codons(genetic_code):
            return "Alternative"

        return "Ambiguous-Definite"

    if status == "possible":
        return "Possible"

    return "None"


# ============================================================
# ORF FINDER
# ============================================================

def resolve_start_codons(start_codons=None, start_codon_policy=None):
    """Return the start-codon set from an explicit set or a named policy."""
    if start_codons is not None:
        return set(start_codons)
    policy = start_codon_policy or DEFAULT_START_CODON_POLICY
    if policy not in START_CODON_POLICIES:
        raise ValueError(
            f"Unknown start codon policy '{policy}'. "
            f"Use one of: {', '.join(START_CODON_POLICIES)}."
        )
    return set(START_CODON_POLICIES[policy])
