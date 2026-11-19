"""The RNA source of the session pipeline (selected ORF or sequence)."""

from __future__ import annotations

from analyses._rna.bases import orf_genetic_code, _orf_value
from analyses._rna.summary import analyze_rna


def _clean_sequence(sequence):
    """Normalize a sequence without changing IUPAC symbols."""

    return (
        str(sequence or "")
        .upper()
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .replace("\t", "")
    )


def _get_pipeline_context(pipeline):
    """Return the pipeline context regardless of pipeline implementation."""

    if pipeline is None:
        return None

    return getattr(pipeline, "context", pipeline)


def _get_pipeline_selected_orf(pipeline):
    """Get the selected ORF from the pipeline if one exists."""

    context = _get_pipeline_context(pipeline)

    if context is None:
        return None

    selected = getattr(context, "selected_orf", None)

    if selected is None:
        selected = getattr(pipeline, "selected_orf", None)

    return selected


def resolve_rna_source(dna, selected_orf=None, pipeline=None):
    """
    Resolve the biological input for RNA analysis.

    Priority:
        1. Explicit selected_orf
        2. selected_orf stored in pipeline
        3. Complete DNA sequence

    This is important because users are allowed to enter RNA Analysis
    directly without selecting an ORF first.
    """

    full_dna = _clean_sequence(dna)

    if selected_orf is None:
        selected_orf = _get_pipeline_selected_orf(pipeline)

    if selected_orf is not None:
        orf_sequence = _clean_sequence(
            _orf_value(selected_orf, "sequence", "")
        )

        if orf_sequence:
            return {
                "dna": orf_sequence,
                "source": "Selected ORF",
                "mode": "orf",
                "orf": selected_orf,
            }

    return {
        "dna": full_dna,
        "source": "Full DNA Sequence",
        "mode": "full_sequence",
        "orf": None,
    }


def _store_rna_in_pipeline(analysis, pipeline):
    """
    Store RNA analysis in whichever pipeline storage API is available.

    This function is deliberately defensive so that RNA works with the
    current pipeline implementation without forcing a specific class design.
    """

    if pipeline is None:
        return

    context = _get_pipeline_context(pipeline)

    if context is None:
        return

    store_method = getattr(
        context,
        "store_rna_analysis",
        None
    )

    if callable(store_method):
        try:
            store_method(
                analysis,
                rna_sequence=analysis.get("rna")
            )
        except TypeError:
            store_method(analysis)

    # Always keep a direct reference when the context allows attributes.
    try:
        context.rna_analysis = analysis
    except Exception:
        pass


def analyze_rna_from_pipeline(
    dna,
    pipeline=None,
    selected_orf=None,
    **analysis_kwargs
):
    """
    Run the original analyze_rna() engine with pipeline-aware input.

    'header' is accepted and ignored for compatibility with main.py.
    It must NOT be forwarded into analyze_rna(), because analyze_rna()
    does not define a header parameter.
    """

    # Compatibility with callers that still pass header=...
    analysis_kwargs.pop("header", None)

    source = resolve_rna_source(
        dna,
        selected_orf=selected_orf,
        pipeline=pipeline,
    )

    if not source["dna"]:
        raise ValueError(
            "No DNA sequence is available for RNA analysis."
        )

    # Audit P1-9: a selected ORF carries the table it was found with, so the
    # codon view must use it unless the caller states another one.
    analysis_kwargs.setdefault(
        "genetic_code",
        orf_genetic_code(source["orf"]),
    )

    analysis = analyze_rna(
        source["dna"],
        **analysis_kwargs
    )

    analysis["pipeline_input"] = {
        "source": source["source"],
        "mode": source["mode"],
        "selected_orf": source["orf"],
        "dna": source["dna"],
        "dna_length": len(source["dna"]),
    }

    _store_rna_in_pipeline(
        analysis,
        pipeline
    )

    return analysis


# ============================================================
# PUBLIC PIPELINE ALIASES
# ============================================================
#
# Keep these aliases so existing main.py / pipeline.py imports
# continue to work. There is only ONE actual UI implementation.
# ============================================================

run_rna_pipeline = analyze_rna_from_pipeline