"""Variants in primer and guide binding sites, drawn after a design (T12, audit P2-4).

The variant dialog records the variants of the VCF contig that the user
mapped onto the loaded sequence (PipelineContext.set_mapped_variants). After
a primer or guide design on that same sequence, this section lists the
variants that fall into the binding sites of the top designs. A design on
another sequence (a manually entered one, an ORF-relative fallback) is not
checked: its positions are not those of the VCF.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, cast

from core.api import binding_features
from core.textformat import columns
from ui import style


def print_variant_check(pipeline: Any, sequence: str, module: str, data: Any,
                        *, top_n: int = 5) -> Optional[List[Dict[str, Any]]]:
    """Draw the variant check for a design, or nothing when no variants apply.

    Returns the findings, or None when no mapped variants exist for this
    sequence.
    """
    from analyses.variant_analysis import check_features_against_variants

    context = getattr(pipeline, "context", None)
    lookup = getattr(context, "variants_for_sequence", None)
    mapping = cast(Optional[Dict[str, Any]], lookup(sequence) if callable(lookup) else None)
    if mapping is None:
        return None
    features = binding_features(module, data, top_n=top_n)
    findings = check_features_against_variants(features, mapping["variants"], chromosome=mapping["contig"])
    noun = "primer pairs" if module == "primer" else "guides"
    style.section("VARIANTS IN BINDING SITES")
    style.field("VCF contig", mapping["contig"])
    style.field("Variants on the contig", len(mapping["variants"]))
    style.field("Designs checked", f"top {top_n} {noun}")
    if not findings:
        style.note("No variant lies in a binding site of these designs.")
        return findings
    widths = [22, 16, 12, 22]
    print(columns(["Design", "Variant", "Severity", "Region"], widths))
    for item in findings:
        print(columns([
            item["feature"],
            f"{item['position']} {item['reference']}>{item['alternative']}",
            item["severity"],
            item["region"],
        ], widths))
    style.warning("A variant in a binding site can reduce binding or cutting; HIGH means "
                  "the primer's 3' end or the guide's seed or PAM.")
    return findings
