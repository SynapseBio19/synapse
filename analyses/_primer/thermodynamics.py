"""Melting temperature and duplex thermodynamics.

Wallace and salt-adjusted screening formulas, the SantaLucia
nearest-neighbor model with its parameter sets, salt and concentration
terms, the backward-compatible ``calculate_tm_nearest_neighbor`` and the
condition presets."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import log, log10, sqrt
from typing import Literal, Optional, overload

from analyses._primer.sequence import count_gc_clamp, normalize_dna_sequence, reverse_complement
from analyses._primer.models import PrimerDesignConfig


# Names of the nearest-neighbor Tm model. ProfessionalPrimerConfig accepts
# "nearest-neighbor", so every dispatch accepts it too (it raised
# "Unknown Tm method" in candidate generation before 2026-11-12).
NEAREST_NEIGHBOR_TM_METHODS = frozenset(
    {"nearest_neighbor", "nearest-neighbor", "nn", "santalucia", "santa_lucia"}
)

# SantaLucia 1998 DNA/DNA nearest-neighbor parameters.
# Values: (delta_H kcal/mol, delta_S cal/(mol*K)).
# Keys are written as 5'-XY-3' / 3'-X'Y'-5' dinucleotide stacks.
NN_PARAMETERS = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

# SantaLucia initiation terms.
NN_INITIATION = (0.2, -5.7)
NN_AT_TERMINAL = (2.2, 6.9)
NN_GC_TERMINAL = (0.0, -0.1)
NN_SYMMETRY_DS = -1.4
GAS_CONSTANT = 1.98720425864083  # cal/(K*mol)


# ============================================================================
# THERMODYNAMICS
# ============================================================================

# Superseded definition of calculate_tm_nearest_neighbor() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


def calculate_tm_wallace(sequence: str) -> float:
    """Wallace-rule screening estimate: Tm = 2(A+T) + 4(G+C)."""
    sequence = normalize_dna_sequence(sequence)
    at = sequence.count("A") + sequence.count("T")
    gc = sequence.count("G") + sequence.count("C")
    return float(2 * at + 4 * gc)


def calculate_tm_salt_adjusted(sequence: str, salt_mM: float = 50.0) -> float:
    """Empirical monovalent-salt-adjusted Tm estimate."""
    sequence = normalize_dna_sequence(sequence)
    if salt_mM <= 0:
        raise ValueError("salt_mM must be greater than zero.")
    n = len(sequence)
    gc = sequence.count("G") + sequence.count("C")
    salt_molar = salt_mM / 1000.0
    return 64.9 + 41.0 * (gc - 16.4) / n + 16.6 * log10(salt_molar)


def calculate_tm(
    sequence: str,
    method: str = "nearest_neighbor",
    salt_mM: float = 50.0,
    mg_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_concentration_nM: float = 250.0,
    include_symmetry_correction: bool = True,
) -> float:
    """Calculate Tm using Wallace, salt-adjusted, or nearest-neighbor model."""
    method = method.lower().strip()
    if method == "wallace":
        return calculate_tm_wallace(sequence)
    if method in {"salt", "salt_adjusted"}:
        return calculate_tm_salt_adjusted(sequence, salt_mM)
    if method in NEAREST_NEIGHBOR_TM_METHODS:
        return float(calculate_tm_nearest_neighbor(
            sequence,
            salt_mM=salt_mM,
            mg_mM=mg_mM,
            dntp_mM=dntp_mM,
            primer_concentration_nM=primer_concentration_nM,
            include_symmetry_correction=include_symmetry_correction,
        ))
    raise ValueError(
        "Unknown Tm method. Use 'nearest_neighbor', 'wallace', or 'salt_adjusted'."
    )


def get_thermodynamic_conditions(config: PrimerDesignConfig) -> dict:
    """Return the configured reaction conditions in a serializable dictionary."""
    return {
        "salt_mM": config.salt_mM,
        "mg_mM": config.mg_mM,
        "dntp_mM": config.dntp_mM,
        "primer_concentration_nM": config.primer_concentration_nM,
        "include_symmetry_correction": config.include_symmetry_correction,
    }
# ============================================================================
# SYNAPSE — PRIMER DESIGN ENGINE V3 THERMODYNAMIC EXTENSION
# ============================================================================
#
# This extension deliberately preserves the complete V2 engine above.
#
# The original candidate generation, scoring, structure heuristics,
# specificity screening, pair ranking, reporting, structured API and
# regression tests remain available.
#
# V3 adds a real nearest-neighbor DNA/DNA thermodynamic screening layer
# with explicit experimental-condition parameters:
#
#     - Na+ / monovalent salt
#     - Mg2+
#     - dNTP concentration
#     - primer concentration
#
# The nearest-neighbor model is based on published DNA/DNA thermodynamic
# parameter sets commonly associated with SantaLucia-style calculations.
#
# IMPORTANT SCIENTIFIC SCOPE
# --------------------------
# This remains a primer-design screening engine, not a replacement for
# experimentally validated commercial thermodynamic software.
#
# The model accounts for duplex thermodynamics and ionic conditions at
# screening level. It does not model every laboratory effect, including
# DMSO, betaine, cosolvents, mismatches, modified bases, non-standard
# buffer systems, polymerase-specific behavior, or sequence-specific
# concentration effects.
#
# The public API remains backward compatible:
#
#     design_primers(...)
#     analyze_primer(...)
#     get_primer_summary(...)
#     get_amplicon(...)
#
# New public functions:
#
#     calculate_tm_nearest_neighbor(...)
#     calculate_duplex_thermodynamics(...)
#     calculate_effective_ionic_strength(...)
#     calculate_terminal_stability(...)
#     calculate_primer_thermodynamics(...)
#     validate_thermodynamic_conditions(...)
#     get_thermodynamic_summary(...)
#
# New configuration:
#
#     ThermodynamicConfig
#
# Existing calls that do not provide a thermodynamic configuration continue
# to use the existing V2 Wallace default.
# ============================================================================



# ============================================================================
# THERMODYNAMIC CONSTANTS
# ============================================================================

GAS_CONSTANT = 1.98720425864083
ABSOLUTE_ZERO_C = 273.15

# SantaLucia-style DNA/DNA nearest-neighbor parameters.
#
# ΔH values are kcal/mol.
# ΔS values are cal/(mol*K).
#
# Keys are written as dinucleotide steps on the strand being evaluated.
# Reverse-complementary steps use the corresponding symmetry-equivalent
# parameterization.
#
# The values are a standard screening set used for DNA duplex Tm work.
NEAREST_NEIGHBOR_PARAMETERS = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}


# Initiation terms.
# These are included explicitly rather than hidden inside the Tm formula.
INITIATION_ENTHALPY = 0.2
INITIATION_ENTROPY = -5.7

# Terminal AT correction.
TERMINAL_AT_ENTHALPY = 2.2
TERMINAL_AT_ENTROPY = 6.9

# Symmetry correction for self-complementary duplexes.
SYMMETRY_ENTROPY = -1.4


# ----------------------------------------------------------------------------
# Consistent nearest-neighbor parameter sets (2026-10-30, Stage 3)
# ----------------------------------------------------------------------------
# Until 2026-10-29 Synapse combined the SantaLucia (1998) stack parameters
# with the SantaLucia & Hicks (2004) initiation terms, i.e. two different
# published parameter sets. The default is now the internally consistent
# SantaLucia (1998) unified set with terminal-specific initiation, the set
# used by Primer3's recommended default (PRIMER_TM_FORMULA=1). The previous
# combination remains available as "synapse_legacy" for reproducibility.
#
# Values: (delta_H kcal/mol, delta_S cal/(mol*K)).
_SANTALUCIA_2004_STACKS = dict(NEAREST_NEIGHBOR_PARAMETERS)
_SANTALUCIA_2004_STACKS.update({
    "AA": (-7.6, -21.3),
    "TT": (-7.6, -21.3),
})

NN_PARAMETER_SETS = {
    "santalucia_1998": {
        "reference": "SantaLucia JR (1998) PNAS 95:1460-1465 (unified parameters)",
        "stacks": dict(NEAREST_NEIGHBOR_PARAMETERS),
        "initiation": (0.0, 0.0),
        "terminal": {"GC": (0.1, -2.8), "AT": (2.3, 4.1)},
        "salt_phosphates": "n_minus_1",
        "self_complementary_concentration_divisor": 1.0,
    },
    "santalucia_2004": {
        "reference": "SantaLucia JR, Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415-440",
        "stacks": _SANTALUCIA_2004_STACKS,
        "initiation": (0.2, -5.7),
        "terminal": {"GC": (0.0, 0.0), "AT": (2.2, 6.9)},
        "salt_phosphates": "n_minus_1",
        "self_complementary_concentration_divisor": 1.0,
    },
    "synapse_legacy": {
        "reference": "Synapse before 2026-10-30: 1998 stacks + 2004 initiation (mixed)",
        "stacks": dict(NEAREST_NEIGHBOR_PARAMETERS),
        "initiation": (INITIATION_ENTHALPY, INITIATION_ENTROPY),
        "terminal": {"GC": (0.0, 0.0), "AT": (TERMINAL_AT_ENTHALPY, TERMINAL_AT_ENTROPY)},
        "salt_phosphates": "n",
        "self_complementary_concentration_divisor": 2.0,
    },
}

DEFAULT_NN_PARAMETER_SET = "santalucia_1998"

# Divalent-cation handling. "monovalent_equivalent" now uses the published
# conversion of von Ahsen et al. (2001), Clin Chem 47:1956-1961, as used by
# Primer3 (oligotm.c divalent_to_monovalent) and Biopython (salt_correction):
#     [Na_eq] (mM) = [Na+] + 120 * sqrt([Mg2+] - [dNTP]),  0 if dNTP >= Mg2+.
# "synapse_legacy_4x" is the previous ad-hoc rule [Na+] + 4 * ([Mg2+] - [dNTP]).
SALT_MODELS = ("monovalent_equivalent", "na_only", "synapse_legacy_4x")


def _nn_parameter_set(name: str) -> dict:
    try:
        return NN_PARAMETER_SETS[name]
    except KeyError:
        raise ValueError(
            f"Unknown nearest-neighbor parameter set '{name}'. "
            f"Use one of: {', '.join(NN_PARAMETER_SETS)}."
        ) from None


def _salt_phosphate_count(sequence: str, parameter_set: str) -> int:
    """Number of phosphates per strand used in the SantaLucia salt term (N - 1)."""
    if _nn_parameter_set(parameter_set)["salt_phosphates"] == "n":
        return len(sequence)
    return len(sequence) - 1


# ============================================================================
# THERMODYNAMIC DATA MODELS
# ============================================================================

@dataclass(frozen=True)
class ThermodynamicConfig:
    """
    Experimental-condition configuration for primer Tm screening.

    Concentrations are expressed in common laboratory units:
        Na+      mM
        Mg2+     mM
        dNTP     mM
        primer   nM

    The defaults are deliberately conservative and are intended as a
    documented screening condition rather than a universal PCR recipe.
    """

    sodium_mM: float = 50.0
    magnesium_mM: float = 1.5
    dntp_mM: float = 0.2
    primer_nM: float = 250.0

    dna_concentration_mode: str = "duplex"
    temperature_reference_C: float = 25.0

    salt_model: str = "monovalent_equivalent"
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET

    def validate(self) -> None:
        validate_thermodynamic_conditions(self)


@dataclass(frozen=True)
class DuplexThermodynamics:
    """
    Thermodynamic result for one DNA/DNA duplex.

    ΔH is kcal/mol.
    ΔS is cal/(mol*K).
    ΔG is kcal/mol at the reported reference temperature.
    """

    sequence: str
    delta_h_kcal_mol: float
    delta_s_cal_mol_k: float
    delta_g_kcal_mol: float
    tm_C: float

    sodium_mM: float
    magnesium_mM: float
    dntp_mM: float
    primer_nM: float

    effective_salt_mM: float
    concentration_term: float

    self_complementary: bool
    terminal_at_count: int


@dataclass(frozen=True)
class PrimerThermodynamicProfile:
    """
    Thermodynamic profile used by candidate scoring.
    """

    method: str
    tm_C: float
    delta_h_kcal_mol: float
    delta_s_cal_mol_k: float
    delta_g_kcal_mol: float

    effective_salt_mM: float

    sodium_mM: float
    magnesium_mM: float
    dntp_mM: float
    primer_nM: float

    concentration_term: float
    terminal_at_count: int
    self_complementary: bool


# ============================================================================
# THERMODYNAMIC VALIDATION
# ============================================================================

def validate_thermodynamic_conditions(
    conditions: ThermodynamicConfig,
) -> None:
    """Validate all thermodynamic input parameters."""

    if conditions.sodium_mM <= 0:
        raise ValueError(
            "sodium_mM must be greater than zero."
        )

    if conditions.magnesium_mM < 0:
        raise ValueError(
            "magnesium_mM cannot be negative."
        )

    if conditions.dntp_mM < 0:
        raise ValueError(
            "dntp_mM cannot be negative."
        )

    if conditions.primer_nM <= 0:
        raise ValueError(
            "primer_nM must be greater than zero."
        )

    if conditions.temperature_reference_C <= -ABSOLUTE_ZERO_C:
        raise ValueError(
            "temperature_reference_C is below absolute zero."
        )

    if conditions.dna_concentration_mode.lower() not in {
        "duplex",
        "strand",
    }:
        raise ValueError(
            "dna_concentration_mode must be 'duplex' or 'strand'."
        )

    if conditions.salt_model.lower() not in SALT_MODELS:
        raise ValueError(
            f"salt_model must be one of: {', '.join(SALT_MODELS)}."
        )

    _nn_parameter_set(conditions.nn_parameter_set)


# ============================================================================
# EFFECTIVE IONIC CONDITIONS
# ============================================================================

def calculate_effective_ionic_strength(
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    *,
    salt_model: str = "monovalent_equivalent",
) -> float:
    """
    Estimate an effective monovalent-equivalent salt concentration.

    The purpose is to provide a transparent screening approximation when
    both Na+ and Mg2+ are supplied.

    dNTPs reduce free Mg2+ availability. A simple conservative correction is
    applied before converting Mg2+ into a monovalent-equivalent term.

    This is intentionally documented as an approximation; a complete PCR
    buffer speciation model is outside the scope of this engine.
    """

    if sodium_mM <= 0:
        raise ValueError(
            "sodium_mM must be greater than zero."
        )

    if magnesium_mM < 0:
        raise ValueError(
            "magnesium_mM cannot be negative."
        )

    if dntp_mM < 0:
        raise ValueError(
            "dntp_mM cannot be negative."
        )

    model = salt_model.lower()

    if model == "na_only":
        return float(sodium_mM)

    if model not in SALT_MODELS:
        raise ValueError(
            "Unknown salt model."
        )

    # Free Mg2+ after complexation by dNTPs (1:1 approximation).
    free_mg = max(
        0.0,
        magnesium_mM - dntp_mM,
    )

    if model == "synapse_legacy_4x":
        # Pre-2026-10-30 Synapse rule, kept only for reproducibility; it has
        # no published basis and underestimates the Mg2+ effect.
        magnesium_equivalent = 4.0 * free_mg
    else:
        # von Ahsen et al. (2001); identical to Primer3 and Biopython.
        magnesium_equivalent = 120.0 * sqrt(free_mg)

    effective = (
        float(sodium_mM)
        + magnesium_equivalent
    )

    return max(
        0.001,
        effective,
    )


def calculate_salt_correction(
    effective_salt_mM: float,
) -> float:
    """
    Calculate the logarithmic monovalent salt correction term.

    This correction is applied to duplex entropy in the nearest-neighbor
    screening model.
    """

    if effective_salt_mM <= 0:
        raise ValueError(
            "effective_salt_mM must be greater than zero."
        )

    salt_M = effective_salt_mM / 1000.0

    return 0.368 * log(
        salt_M
    )


# ============================================================================
# SEQUENCE-LEVEL THERMODYNAMIC HELPERS
# ============================================================================

def _is_self_complementary(
    sequence: str,
) -> bool:
    """Return True when a sequence equals its reverse complement."""

    sequence = normalize_dna_sequence(
        sequence
    )

    return (
        sequence
        == reverse_complement(sequence)
    )


def _terminal_at_count(
    sequence: str,
) -> int:
    """Count A/T residues at the two duplex termini."""

    sequence = normalize_dna_sequence(
        sequence
    )

    count = 0

    if sequence[0] in "AT":
        count += 1

    if sequence[-1] in "AT":
        count += 1

    return count


def _nearest_neighbor_sum(
    sequence: str,
    parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> tuple[float, float]:
    """
    Sum nearest-neighbor ΔH and ΔS values (1 M NaCl) for one parameter set,
    including initiation, terminal and symmetry terms.
    """
    # A primer's sum is needed for its Tm and again for its duplex profile;
    # it depends only on the two strings, so it is cached (P2-6).
    if isinstance(sequence, str) and isinstance(parameter_set, str):
        return _nearest_neighbor_sum_cached(sequence, parameter_set)
    return _nearest_neighbor_sum_uncached(sequence, parameter_set)


def _nearest_neighbor_sum_uncached(
    sequence: str,
    parameter_set: str,
) -> tuple[float, float]:
    sequence = normalize_dna_sequence(
        sequence
    )

    if len(sequence) < 2:
        raise ValueError(
            "Nearest-neighbor thermodynamics requires at least 2 bases."
        )

    parameters = _nn_parameter_set(parameter_set)
    stacks = parameters["stacks"]

    delta_h, delta_s = parameters["initiation"]

    for left, right in zip(
        sequence,
        sequence[1:],
    ):
        key = left + right

        if key not in stacks:
            raise ValueError(
                f"No nearest-neighbor parameter for dinucleotide {key}."
            )

        step_h, step_s = stacks[key]

        delta_h += step_h
        delta_s += step_s

    for terminal_base in (sequence[0], sequence[-1]):
        kind = "AT" if terminal_base in "AT" else "GC"
        term_h, term_s = parameters["terminal"][kind]
        delta_h += term_h
        delta_s += term_s

    if _is_self_complementary(
        sequence
    ):
        delta_s += SYMMETRY_ENTROPY

    return (
        delta_h,
        delta_s,
    )


_nearest_neighbor_sum_cached = lru_cache(maxsize=65536)(_nearest_neighbor_sum_uncached)


def calculate_duplex_concentration_term(
    primer_nM: float,
    *,
    self_complementary: bool = False,
    mode: str = "duplex",
    parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> float:
    """
    Calculate the concentration-dependent logarithmic entropy term.

    The distinction between duplex and strand concentration conventions is
    retained explicitly so the calculation is inspectable.
    """

    if primer_nM <= 0:
        raise ValueError(
            "primer_nM must be greater than zero."
        )

    mode = mode.lower()

    concentration_M = (
        primer_nM
        * 1e-9
    )

    if mode == "duplex":
        concentration_term = log(
            concentration_M / 4.0
        )
    elif mode == "strand":
        concentration_term = log(
            concentration_M
            / 2.0
        )
    else:
        raise ValueError(
            "mode must be 'duplex' or 'strand'."
        )

    if self_complementary:
        # SantaLucia (1998): CT for self-complementary duplexes (Primer3
        # "Equation A"); the legacy set used CT/2.
        divisor = _nn_parameter_set(parameter_set)[
            "self_complementary_concentration_divisor"
        ]
        concentration_term = log(
            concentration_M / divisor
        )

    return concentration_term


# ============================================================================
# NEAREST-NEIGHBOR Tm
# ============================================================================

# This is the V3 model itself. The public entry point is the
# backward-compatible adapter calculate_tm_nearest_neighbor() further down,
# which accepts the older parameter names and delegates here.
#
# Audit D16: this implementation used to carry the public name and was then
# rebound to the private name and shadowed by the adapter, so the definition a
# reader found first was not the one that ran. Naming it privately here leaves
# exactly one definition per name without changing which code executes.
def _calculate_tm_nearest_neighbor_raw(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> float:
    """
    Calculate a nearest-neighbor DNA/DNA duplex Tm.

    Parameters:
        sodium_mM:
            Monovalent sodium concentration in mM.

        magnesium_mM:
            Mg2+ concentration in mM.

        dntp_mM:
            Total dNTP concentration in mM.

        primer_nM:
            Primer concentration in nM.

        dna_concentration_mode:
            'duplex' or 'strand'.

        salt_model:
            'monovalent_equivalent' or 'na_only'.

    Returns:
        Melting temperature in degrees Celsius.

    Scientific note:
        This is a sequence-level thermodynamic screening model. It is not
        intended to replace experimental calibration or a complete PCR
        buffer model.
    """

    conditions = ThermodynamicConfig(
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        dna_concentration_mode=dna_concentration_mode,
        salt_model=salt_model,
        nn_parameter_set=nn_parameter_set,
    )

    conditions.validate()

    sequence = normalize_dna_sequence(
        sequence
    )

    delta_h, delta_s = (
        _nearest_neighbor_sum(
            sequence,
            nn_parameter_set,
        )
    )

    effective_salt = (
        calculate_effective_ionic_strength(
            sodium_mM=sodium_mM,
            magnesium_mM=magnesium_mM,
            dntp_mM=dntp_mM,
            salt_model=salt_model,
        )
    )

    salt_correction = (
        calculate_salt_correction(
            effective_salt
        )
    )

    corrected_delta_s = (
        delta_s
        + _salt_phosphate_count(sequence, nn_parameter_set)
        * salt_correction
    )

    self_complementary = (
        _is_self_complementary(
            sequence
        )
    )

    concentration_term = (
        calculate_duplex_concentration_term(
            primer_nM,
            self_complementary=self_complementary,
            mode=dna_concentration_mode,
            parameter_set=nn_parameter_set,
        )
    )

    denominator = (
        corrected_delta_s
        + GAS_CONSTANT
        * concentration_term
    )

    if denominator == 0:
        raise ZeroDivisionError(
            "Thermodynamic denominator is zero."
        )

    tm_K = (
        (delta_h * 1000.0)
        / denominator
    )

    tm_C = (
        tm_K
        - ABSOLUTE_ZERO_C
    )

    return float(
        tm_C
    )


# ============================================================================
# FULL DUPLEX THERMODYNAMICS
# ============================================================================

def calculate_duplex_thermodynamics(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    reference_temperature_C: float = 25.0,
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> DuplexThermodynamics:
    """
    Return the full nearest-neighbor thermodynamic profile.
    """
    arguments = (
        sequence, sodium_mM, magnesium_mM, dntp_mM, primer_nM,
        dna_concentration_mode, salt_model, reference_temperature_C,
        nn_parameter_set,
    )
    try:
        hash(arguments)
    except TypeError:
        return _duplex_thermodynamics(*arguments)
    return _duplex_thermodynamics_cached(*arguments)


def _duplex_thermodynamics(
    sequence, sodium_mM, magnesium_mM, dntp_mM, primer_nM,
    dna_concentration_mode, salt_model, reference_temperature_C,
    nn_parameter_set,
) -> DuplexThermodynamics:
    conditions = ThermodynamicConfig(
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        dna_concentration_mode=dna_concentration_mode,
        temperature_reference_C=reference_temperature_C,
        salt_model=salt_model,
        nn_parameter_set=nn_parameter_set,
    )

    conditions.validate()

    sequence = normalize_dna_sequence(
        sequence
    )

    delta_h, delta_s = (
        _nearest_neighbor_sum(
            sequence,
            nn_parameter_set,
        )
    )

    effective_salt = (
        calculate_effective_ionic_strength(
            sodium_mM=sodium_mM,
            magnesium_mM=magnesium_mM,
            dntp_mM=dntp_mM,
            salt_model=salt_model,
        )
    )

    salt_correction = (
        calculate_salt_correction(
            effective_salt
        )
    )

    corrected_delta_s = (
        delta_s
        + _salt_phosphate_count(sequence, nn_parameter_set)
        * salt_correction
    )

    self_complementary = (
        _is_self_complementary(
            sequence
        )
    )

    concentration_term = (
        calculate_duplex_concentration_term(
            primer_nM,
            self_complementary=self_complementary,
            mode=dna_concentration_mode,
            parameter_set=nn_parameter_set,
        )
    )

    denominator = (
        corrected_delta_s
        + GAS_CONSTANT
        * concentration_term
    )

    tm_K = (
        (delta_h * 1000.0)
        / denominator
    )

    tm_C = (
        tm_K
        - ABSOLUTE_ZERO_C
    )

    reference_K = (
        reference_temperature_C
        + ABSOLUTE_ZERO_C
    )

    delta_g = (
        delta_h
        - (
            reference_K
            * corrected_delta_s
            / 1000.0
        )
    )

    return DuplexThermodynamics(
        sequence=sequence,
        delta_h_kcal_mol=round(
            delta_h,
            4,
        ),
        delta_s_cal_mol_k=round(
            corrected_delta_s,
            4,
        ),
        delta_g_kcal_mol=round(
            delta_g,
            4,
        ),
        tm_C=round(
            tm_C,
            4,
        ),
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        effective_salt_mM=round(
            effective_salt,
            4,
        ),
        concentration_term=round(
            concentration_term,
            6,
        ),
        self_complementary=self_complementary,
        terminal_at_count=_terminal_at_count(
            sequence
        ),
    )


# ============================================================================
# TERMINAL STABILITY
# ============================================================================

def calculate_terminal_stability(
    sequence: str,
    window: int = 5,
) -> dict:
    """
    Analyze the thermodynamic composition of the primer 3' terminus.

    This is complementary to the existing GC-clamp calculation. The result
    describes the terminal window rather than merely counting G/C residues.
    """

    sequence = normalize_dna_sequence(
        sequence
    )

    if window <= 0:
        raise ValueError(
            "window must be greater than zero."
        )

    terminal = sequence[
        -min(window, len(sequence)):
    ]

    gc = (
        terminal.count("G")
        + terminal.count("C")
    )

    at = (
        terminal.count("A")
        + terminal.count("T")
    )

    return {
        "window": len(terminal),
        "sequence": terminal,
        "gc_count": gc,
        "at_count": at,
        "gc_percent": (
            round(
                100.0 * gc / len(terminal),
                2,
            )
            if terminal
            else 0.0
        ),
        "gc_clamp": count_gc_clamp(
            sequence
        ),
    }


# ============================================================================
# PRIMER THERMODYNAMIC PROFILE
# ============================================================================

def calculate_primer_thermodynamics(
    sequence: str,
    *,
    conditions: ThermodynamicConfig = ThermodynamicConfig(),
) -> PrimerThermodynamicProfile:
    """
    Calculate the complete thermodynamic profile for one primer.
    """

    conditions.validate()

    duplex = (
        calculate_duplex_thermodynamics(
            sequence,
            sodium_mM=conditions.sodium_mM,
            magnesium_mM=conditions.magnesium_mM,
            dntp_mM=conditions.dntp_mM,
            primer_nM=conditions.primer_nM,
            dna_concentration_mode=(
                conditions.dna_concentration_mode
            ),
            salt_model=conditions.salt_model,
            reference_temperature_C=(
                conditions.temperature_reference_C
            ),
            nn_parameter_set=conditions.nn_parameter_set,
        )
    )

    return PrimerThermodynamicProfile(
        method="nearest_neighbor",
        tm_C=duplex.tm_C,
        delta_h_kcal_mol=(
            duplex.delta_h_kcal_mol
        ),
        delta_s_cal_mol_k=(
            duplex.delta_s_cal_mol_k
        ),
        delta_g_kcal_mol=(
            duplex.delta_g_kcal_mol
        ),
        effective_salt_mM=(
            duplex.effective_salt_mM
        ),
        sodium_mM=(
            conditions.sodium_mM
        ),
        magnesium_mM=(
            conditions.magnesium_mM
        ),
        dntp_mM=(
            conditions.dntp_mM
        ),
        primer_nM=(
            conditions.primer_nM
        ),
        concentration_term=(
            duplex.concentration_term
        ),
        terminal_at_count=(
            duplex.terminal_at_count
        ),
        self_complementary=(
            duplex.self_complementary
        ),
    )


# Primer design asks for the same duplex profile twice per candidate (once
# for the candidate's Tm, once for its thermodynamic profile). The result is
# an immutable value that depends only on the arguments, so it is cached;
# typed=True keeps 50 and 50.0 apart because the record stores the value as
# given (P2-6). An invalid input raises every time and is never cached.
_duplex_thermodynamics_cached = lru_cache(maxsize=65536, typed=True)(_duplex_thermodynamics)


# ============================================================================
# THERMODYNAMIC METHOD DISPATCH
# ============================================================================

def calculate_tm_advanced(
    sequence: str,
    *,
    method: str = "nearest_neighbor",
    thermodynamic_config: ThermodynamicConfig = (
        ThermodynamicConfig()
    ),
    salt_mM: float = 50.0,
) -> float:
    """
    Unified Tm dispatcher for the professional engine.
    """

    method = method.lower().strip()

    if method == "wallace":
        return calculate_tm_wallace(
            sequence
        )

    if method in {
        "salt",
        "salt_adjusted",
    }:
        return calculate_tm_salt_adjusted(
            sequence,
            salt_mM,
        )

    if method in {
        "nearest_neighbor",
        "nearest-neighbor",
        "nn",
        "santalucia",
    }:
        return calculate_tm_nearest_neighbor(
            sequence,
            sodium_mM=(
                thermodynamic_config.sodium_mM
            ),
            magnesium_mM=(
                thermodynamic_config.magnesium_mM
            ),
            dntp_mM=(
                thermodynamic_config.dntp_mM
            ),
            primer_nM=(
                thermodynamic_config.primer_nM
            ),
            dna_concentration_mode=(
                thermodynamic_config.dna_concentration_mode
            ),
            salt_model=(
                thermodynamic_config.salt_model
            ),
            nn_parameter_set=(
                thermodynamic_config.nn_parameter_set
            ),
        )

    raise ValueError(
        "Unknown Tm method. "
        "Use 'wallace', 'salt_adjusted', "
        "or 'nearest_neighbor'."
    )


# ============================================================================
# THERMODYNAMIC COMPARISON UTILITIES
# ============================================================================

def compare_tm_methods(
    sequence: str,
    *,
    conditions: ThermodynamicConfig = (
        ThermodynamicConfig()
    ),
) -> dict:
    """
    Compare Wallace, salt-adjusted and nearest-neighbor Tm estimates.
    """

    sequence = normalize_dna_sequence(
        sequence
    )

    wallace = (
        calculate_tm_wallace(
            sequence
        )
    )

    salt_adjusted = (
        calculate_tm_salt_adjusted(
            sequence,
            conditions.sodium_mM,
        )
    )

    nearest_neighbor = (
        calculate_tm_nearest_neighbor(
            sequence,
            sodium_mM=conditions.sodium_mM,
            magnesium_mM=conditions.magnesium_mM,
            dntp_mM=conditions.dntp_mM,
            primer_nM=conditions.primer_nM,
            dna_concentration_mode=(
                conditions.dna_concentration_mode
            ),
            salt_model=conditions.salt_model,
        )
    )

    return {
        "sequence": sequence,
        "wallace_tm_C": round(
            wallace,
            2,
        ),
        "salt_adjusted_tm_C": round(
            salt_adjusted,
            2,
        ),
        "nearest_neighbor_tm_C": round(
            nearest_neighbor,
            2,
        ),
        "wallace_minus_nn_C": round(
            wallace - nearest_neighbor,
            2,
        ),
        "salt_minus_nn_C": round(
            salt_adjusted - nearest_neighbor,
            2,
        ),
    }


# ============================================================================
# DOCUMENTED EXAMPLE CONFIGURATIONS
# ============================================================================

DEFAULT_PROFESSIONAL_THERMODYNAMIC_CONFIG = (
    ThermodynamicConfig(
        sodium_mM=50.0,
        magnesium_mM=1.5,
        dntp_mM=0.2,
        primer_nM=250.0,
        dna_concentration_mode="duplex",
        salt_model="monovalent_equivalent",
    )
)


HIGH_SALT_SCREENING_CONFIG = (
    ThermodynamicConfig(
        sodium_mM=100.0,
        magnesium_mM=1.5,
        dntp_mM=0.2,
        primer_nM=250.0,
        dna_concentration_mode="duplex",
        salt_model="monovalent_equivalent",
    )
)


LOW_SALT_SCREENING_CONFIG = (
    ThermodynamicConfig(
        sodium_mM=20.0,
        magnesium_mM=1.0,
        dntp_mM=0.2,
        primer_nM=250.0,
        dna_concentration_mode="duplex",
        salt_model="monovalent_equivalent",
    )
)


# ============================================================================
# EXPERIMENTAL-CONDITION SUMMARY
# ============================================================================

def get_thermodynamic_summary(
    conditions: ThermodynamicConfig,
) -> dict:
    """
    Return a structured summary of all thermodynamic assumptions.
    """

    conditions.validate()

    effective_salt = (
        calculate_effective_ionic_strength(
            sodium_mM=conditions.sodium_mM,
            magnesium_mM=conditions.magnesium_mM,
            dntp_mM=conditions.dntp_mM,
            salt_model=conditions.salt_model,
        )
    )

    return {
        "model": (
            "DNA/DNA nearest-neighbor"
        ),
        "sodium_mM": (
            conditions.sodium_mM
        ),
        "magnesium_mM": (
            conditions.magnesium_mM
        ),
        "dntp_mM": (
            conditions.dntp_mM
        ),
        "primer_nM": (
            conditions.primer_nM
        ),
        "effective_salt_mM": (
            round(
                effective_salt,
                4,
            )
        ),
        "salt_model": (
            conditions.salt_model
        ),
        "dna_concentration_mode": (
            conditions.dna_concentration_mode
        ),
        "reference_temperature_C": (
            conditions.temperature_reference_C
        ),
        "scientific_scope": (
            "Sequence-level thermodynamic screening; "
            "not a complete PCR buffer-speciation model."
        ),
    }
# ============================================================================
# BACKWARD-COMPATIBILITY ADAPTER
# ============================================================================
#
# The original Synapse primer engine already exposed a nearest-neighbor API
# with legacy keyword names. The adapter below preserves those names while
# keeping the new V3 thermodynamic implementation underneath.
# ============================================================================

@dataclass(frozen=True)
class NearestNeighborTmResult:
    """Backward-compatible detailed nearest-neighbor Tm result."""

    tm: float
    model: str
    sequence: str
    delta_h_kcal_mol: float
    delta_s_cal_mol_k: float
    effective_salt_mM: float
    sodium_mM: float
    magnesium_mM: float
    dntp_mM: float
    primer_concentration_nM: float
    symmetry_correction_applied: bool = False
    explanation: str = (
        "Nearest-neighbor DNA/DNA thermodynamic screening model."
    )

    @property
    def delta_h_kcal(self) -> float:
        return self.delta_h_kcal_mol

    @property
    def delta_s_cal(self) -> float:
        return self.delta_s_cal_mol_k

    @property
    def concentration_nM(self) -> float:
        return self.primer_concentration_nM

    @property
    def salt_mM(self) -> float:
        return self.sodium_mM

    @property
    def mg_mM(self) -> float:
        return self.magnesium_mM


# The return type follows return_details (typing only; one implementation).
@overload
def calculate_tm_nearest_neighbor(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    salt_mM: Optional[float] = None,
    mg_mM: Optional[float] = None,
    primer_concentration_nM: Optional[float] = None,
    include_symmetry_correction: bool = True,
    return_details: Literal[False] = False,
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> float: ...


@overload
def calculate_tm_nearest_neighbor(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    salt_mM: Optional[float] = None,
    mg_mM: Optional[float] = None,
    primer_concentration_nM: Optional[float] = None,
    include_symmetry_correction: bool = True,
    return_details: Literal[True],
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> NearestNeighborTmResult: ...


@overload
def calculate_tm_nearest_neighbor(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    salt_mM: Optional[float] = None,
    mg_mM: Optional[float] = None,
    primer_concentration_nM: Optional[float] = None,
    include_symmetry_correction: bool = True,
    return_details: bool = False,
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> float | NearestNeighborTmResult: ...


def calculate_tm_nearest_neighbor(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    salt_mM: Optional[float] = None,
    mg_mM: Optional[float] = None,
    primer_concentration_nM: Optional[float] = None,
    include_symmetry_correction: bool = True,
    return_details: bool = False,
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> float | NearestNeighborTmResult:
    """
    Backward-compatible nearest-neighbor Tm API.

    Supported parameter names include both the V3 names and the names used
    by earlier Synapse primer-engine tests:

        sodium_mM
        magnesium_mM / mg_mM
        dntp_mM
        primer_nM / primer_concentration_nM
        salt_mM
        return_details
    """

    if salt_mM is not None:
        sodium_mM = float(salt_mM)

    if mg_mM is not None:
        magnesium_mM = float(mg_mM)

    if primer_concentration_nM is not None:
        primer_nM = float(
            primer_concentration_nM
        )

    if not include_symmetry_correction:
        # The symmetry term only matters for self-complementary oligos; it
        # cannot be switched off in the consistent parameter sets.
        import logging
        logging.getLogger("synapse.primer").info(
            "include_symmetry_correction=False is ignored; the symmetry term "
            "is part of the nearest-neighbor model."
        )

    tm = _calculate_tm_nearest_neighbor_raw(
        sequence,
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        dna_concentration_mode=dna_concentration_mode,
        salt_model=salt_model,
        nn_parameter_set=nn_parameter_set,
    )

    if not return_details:
        return tm

    duplex = calculate_duplex_thermodynamics(
        sequence,
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        dna_concentration_mode=dna_concentration_mode,
        salt_model=salt_model,
        nn_parameter_set=nn_parameter_set,
    )

    return NearestNeighborTmResult(
        tm=duplex.tm_C,
        model="nearest_neighbor",
        sequence=duplex.sequence,
        delta_h_kcal_mol=(
            duplex.delta_h_kcal_mol
        ),
        delta_s_cal_mol_k=(
            duplex.delta_s_cal_mol_k
        ),
        effective_salt_mM=(
            duplex.effective_salt_mM
        ),
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_concentration_nM=primer_nM,
        symmetry_correction_applied=(
            duplex.self_complementary
            and include_symmetry_correction
        ),
        explanation=(
            "Nearest-neighbor DNA/DNA thermodynamic screening "
            "with ionic and concentration corrections."
        ),
    )
