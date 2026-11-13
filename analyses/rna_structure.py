"""Optional energy-based RNA secondary structure (ViennaRNA RNAfold).

The RNA module's hairpin screen (``analyses.rna.find_rna_hairpins``) only
looks for complementary stems around a loop; it computes no free energy.
When the ViennaRNA Python package is installed, this module adds what the
screen cannot provide:

- the minimum free energy (MFE) structure of the whole sequence and its
  energy, computed with ViennaRNA's default model (Turner 2004
  nearest-neighbor parameters, 37 degrees C unless another temperature is
  given; Lorenz et al. 2011, Algorithms Mol Biol 6:26);
- the hairpins of that MFE structure;
- for each screen candidate, the free energy of its stem-loop evaluated
  on its own, and how many of its pairs the MFE structure contains.

ViennaRNA is optional, like Biopython and primer3-py: nothing here is
imported at module load, and every function reports ``available: False``
with the reason when the package cannot be loaded. The package is not
redistributed with Synapse (its licence permits research, educational and
commercial use without redistribution for a fee; see
docs/TOOL_AND_DATA_LEDGER.md).

An MFE structure is a thermodynamic prediction for one molecule in
isolation, not an observed structure; suboptimal structures, proteins,
ions and kinetics are not considered.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

ENERGY_MODEL = "Turner 2004 nearest-neighbor parameters (ViennaRNA default)"
REFERENCE = "Lorenz R. et al. (2011) ViennaRNA Package 2.0. Algorithms Mol Biol 6:26"
DEFAULT_MAX_LENGTH = 3000
"""Longest sequence folded by default. MFE folding takes cubic time: on the
development machine 2,000 nt took about 2 s and 4,000 nt about 14 s; callers
may raise the limit explicitly."""

_STANDARD = frozenset("ACGU")
_PAIRS = ({"A", "U"}, {"G", "C"}, {"G", "U"})


def viennarna_status() -> Dict[str, Any]:
    """Whether the ViennaRNA package can be loaded, and its version."""
    try:
        import RNA  # type: ignore[import-not-found]
    except Exception as error:  # ImportError, or OSError when a DLL is blocked
        return {"available": False, "version": None,
                "reason": f"ViennaRNA is not available ({type(error).__name__}: {error})."}
    return {"available": True, "version": getattr(RNA, "__version__", "unknown"), "reason": None}


def _pair_table(structure: str) -> Dict[int, int]:
    stack: List[int] = []
    table: Dict[int, int] = {}
    for index, symbol in enumerate(structure):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")":
            partner = stack.pop()
            table[partner] = index
            table[index] = partner
    return table


def mfe_hairpins(structure: str) -> List[Dict[str, Any]]:
    """Hairpins of a dot-bracket structure, 1-based and inclusive.

    A hairpin is a pair (i, j) with no pair between them; its stem is the
    run of directly stacked pairs closing it (bulges and interior loops end
    the stem). ``start``/``end`` are the outermost stem pair.
    """
    table = _pair_table(structure)
    hairpins = []
    for i, j in sorted(table.items()):
        if i >= j or any(k in table for k in range(i + 1, j)):
            continue
        outer_i, outer_j = i, j
        while outer_i - 1 >= 0 and table.get(outer_i - 1) == outer_j + 1:
            outer_i -= 1
            outer_j += 1
        hairpins.append({
            "start": outer_i + 1,
            "end": outer_j + 1,
            "stem_length": i - outer_i + 1,
            "loop_start": i + 2,
            "loop_end": j,
            "loop_length": j - i - 1,
        })
    return hairpins


def _candidate_structure(candidate: Dict[str, Any]) -> str:
    """Dot-bracket of a screen candidate's stem-loop; non-pairing positions stay unpaired."""
    left, right = candidate["left_stem"], candidate["right_stem"]
    length = candidate["stem_length"]
    left_part = "".join(
        "(" if {left[k], right[length - 1 - k]} in _PAIRS else "." for k in range(length))
    right_part = "".join(")" if symbol == "(" else "." for symbol in reversed(left_part))
    return left_part + "." * candidate["loop_length"] + right_part


def annotate_hairpin_candidates(rna: str, candidates: Sequence[Dict[str, Any]],
                                structure: str, *, temperature: float = 37.0) -> List[Dict[str, Any]]:
    """Energy and MFE support of the hairpin screen's candidates.

    ``isolated_delta_g`` is the free energy (kcal/mol) of the candidate's
    stem-loop evaluated on its own subsequence; a positive value means the
    stem-loop is not favorable even in isolation. ``pairs_in_mfe`` counts
    the candidate's base pairs that the MFE structure of the whole sequence
    contains.
    """
    import RNA  # type: ignore[import-not-found]

    model = RNA.md()
    model.temperature = temperature
    table = _pair_table(structure)
    annotated = []
    for candidate in candidates:
        start0 = candidate["start"] - 1
        sub = rna[start0:candidate["end"]]
        stem_structure = _candidate_structure(candidate)
        energy = RNA.fold_compound(sub, model).eval_structure(stem_structure)
        right_start = start0 + candidate["stem_length"] + candidate["loop_length"]
        pairs = [(start0 + k, right_start + candidate["stem_length"] - 1 - k)
                 for k in range(candidate["stem_length"])
                 if stem_structure[k] == "("]
        annotated.append({
            "start": candidate["start"],
            "end": candidate["end"],
            "stem_length": candidate["stem_length"],
            "loop_length": candidate["loop_length"],
            "isolated_delta_g": round(energy, 2),
            "pairs": len(pairs),
            "pairs_in_mfe": sum(1 for a, b in pairs if table.get(a) == b),
        })
    return annotated


def predict_secondary_structure(rna: str, *, temperature: float = 37.0,
                                max_length: Optional[int] = DEFAULT_MAX_LENGTH,
                                hairpin_candidates: Optional[Sequence[Dict[str, Any]]] = None
                                ) -> Dict[str, Any]:
    """MFE secondary structure of an RNA sequence with ViennaRNA, if installed.

    ``rna`` may be written with T or U; positions are 1-based. Letters other
    than A, C, G and U cannot pair in the energy model; their number is
    reported as a warning. Returns ``available: False`` with the reason when
    ViennaRNA cannot be loaded or the sequence is longer than ``max_length``.
    """
    sequence = "".join(rna.split()).upper().replace("T", "U")
    status = viennarna_status()
    result: Dict[str, Any] = {
        "available": False,
        "tool": "ViennaRNA",
        "version": status["version"],
        "energy_model": ENERGY_MODEL,
        "reference": REFERENCE,
        "temperature_c": temperature,
        "length": len(sequence),
        "warnings": [],
        "limitations": [
            "Thermodynamic prediction for one molecule in isolation (minimum free energy); "
            "suboptimal structures, pseudoknots, proteins, ions and folding kinetics are not considered.",
            "Not an experimentally determined structure.",
        ],
    }
    if not status["available"]:
        result["reason"] = status["reason"]
        return result
    if not sequence:
        result["reason"] = "The sequence is empty."
        return result
    if max_length is not None and len(sequence) > max_length:
        result["reason"] = (f"The sequence has {len(sequence):,} nt; MFE folding is limited to "
                            f"{max_length:,} nt by default because its run time grows with the cube "
                            f"of the length.")
        return result

    import RNA  # type: ignore[import-not-found]

    nonstandard = sum(1 for base in sequence if base not in _STANDARD)
    if nonstandard:
        result["warnings"].append(
            f"{nonstandard} position(s) are not A, C, G or U; they cannot pair in the energy model.")
    model = RNA.md()
    model.temperature = temperature
    structure, energy = RNA.fold_compound(sequence, model).mfe()
    result.update({
        "available": True,
        "reason": None,
        "structure": structure,
        "mfe_kcal_per_mol": round(energy, 2),
        "paired_fraction": round(structure.count("(") * 2 / len(sequence), 4),
        "mfe_hairpins": mfe_hairpins(structure),
    })
    if hairpin_candidates is not None:
        annotated = annotate_hairpin_candidates(sequence, hairpin_candidates, structure,
                                                temperature=temperature)
        result["hairpin_candidates"] = annotated
        result["candidates_favorable"] = sum(1 for item in annotated if item["isolated_delta_g"] < 0)
        result["candidates_in_mfe"] = sum(
            1 for item in annotated if item["pairs"] and item["pairs_in_mfe"] * 2 >= item["pairs"])
    return result


__all__ = [
    "DEFAULT_MAX_LENGTH", "ENERGY_MODEL", "REFERENCE", "annotate_hairpin_candidates",
    "mfe_hairpins", "predict_secondary_structure", "viennarna_status",
]
