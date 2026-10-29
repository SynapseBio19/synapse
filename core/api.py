"""Non-interactive entry points for the Synapse analyses.

Every function here takes explicit parameters, never calls ``input()`` or
``print()`` for its result, and returns an :class:`AnalysisResult` that
records the parameters, the input hash, warnings and limitations. This is
the layer a future graphical or web interface (or a script) should call;
the terminal menus in ``main.py`` remain unchanged.

The functions are thin wrappers: the scientific work is done by the
existing modules in ``analyses/``. Parameter defaults are those of the
underlying modules, so results are identical to the menus' results for the
same settings.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

from core.models import AnalysisResult, SequenceRecord, sha256_text
from core.seqtools import normalize_dna


def _result(module: str, data: Any, sequence: str, description: str, *,
            parameters: Optional[Dict[str, Any]] = None,
            warnings: Optional[List[str]] = None,
            limitations: Optional[List[str]] = None) -> AnalysisResult:
    return AnalysisResult(
        module=module,
        data=data,
        parameters=dict(parameters or {}),
        input_sequence_sha256=sha256_text(sequence) if sequence else None,
        input_description=description,
        warnings=list(warnings or []),
        limitations=list(limitations or []),
    )


def _sequence_of(record_or_sequence: Any) -> tuple[str, str, str]:
    """Return (sequence, description, topology) for a record or a plain string."""
    if isinstance(record_or_sequence, SequenceRecord):
        return (record_or_sequence.sequence, record_or_sequence.identifier,
                record_or_sequence.topology)
    return normalize_dna(str(record_or_sequence)), "sequence", "linear"


# ----------------------------------------------------------------------
# Input
# ----------------------------------------------------------------------
def load_records(path: str) -> List[SequenceRecord]:
    """Load FASTA or GenBank records with provenance."""
    from fileio.loader import load_sequences

    return [SequenceRecord.from_legacy_dict(record, index)
            for index, record in enumerate(load_sequences(path), start=1)]


# ----------------------------------------------------------------------
# Analyses
# ----------------------------------------------------------------------
def sequence_overview(record_or_sequence: Any, *, profile_window: int = 100,
                      profile_step: int = 1) -> AnalysisResult:
    from analyses.sequence import analyze_sequence

    sequence, description, _ = _sequence_of(record_or_sequence)
    data = analyze_sequence(sequence, profile_window=profile_window, profile_step=profile_step)
    return _result("sequence", data, sequence, description,
                   parameters={"profile_window": profile_window, "profile_step": profile_step})


def find_orfs(record_or_sequence: Any, *, min_length: Optional[int] = None,
              start_codon_policy: Optional[str] = None, nested: Optional[bool] = None,
              ranking: Optional[str] = None, include_partial: bool = True,
              genetic_code: int = 1, top_n: int = 5) -> AnalysisResult:
    """ORF search; ``None`` means the module default (see docs/METHODS.md)."""
    from analyses import orf as O

    sequence, description, topology = _sequence_of(record_or_sequence)
    data = O.analyze_orfs_enhanced(
        sequence,
        include_partial=include_partial,
        min_length=O.DEFAULT_MIN_ORF_LENGTH if min_length is None else min_length,
        top_n=top_n,
        start_codon_policy=start_codon_policy,
        nested=O.DEFAULT_REPORT_NESTED if nested is None else nested,
        ranking=ranking or O.DEFAULT_RANKING,
        genetic_code=genetic_code,
    )
    limitations = ["ORFs are sequence features, not evidence of protein-coding genes."]
    if topology == "circular":
        limitations.append("ORFs spanning the origin of a circular sequence are not reported.")
    return _result("orf", data, sequence, description,
                   parameters=data.get("parameters", {}), limitations=limitations)


def analyze_rna(record_or_sequence: Any, *, profile_window: int = 100,
                profile_step: int = 1) -> AnalysisResult:
    from analyses.rna import analyze_rna as _analyze_rna

    sequence, description, _ = _sequence_of(record_or_sequence)
    data = _analyze_rna(sequence, profile_window=profile_window, profile_step=profile_step)
    return _result("rna", data, sequence, description,
                   parameters={"profile_window": profile_window, "profile_step": profile_step})


def analyze_protein_from_cds(cds: str, *, initiator_start: bool = True,
                             genetic_code: int = 1, description: str = "CDS") -> AnalysisResult:
    """Translate a coding sequence (frame +1) and analyse the protein."""
    from analyses.protein import analyze_complete_protein

    sequence = normalize_dna(cds)
    data = analyze_complete_protein("", rna=sequence.replace("T", "U"),
                                    initiator_start=initiator_start, genetic_code=genetic_code)
    return _result("protein", data, sequence, description,
                   parameters={"initiator_start": initiator_start, "genetic_code": genetic_code},
                   limitations=["Physicochemical values are sequence-based estimates."])


def restriction_digest(record_or_sequence: Any, enzymes: Sequence[str], *,
                       topology: Optional[str] = None,
                       include_possible: bool = False) -> AnalysisResult:
    from analyses.restriction import perform_digest

    sequence, description, record_topology = _sequence_of(record_or_sequence)
    topology = topology or record_topology
    data = perform_digest(sequence, list(enzymes), topology=topology,
                          include_possible=include_possible)
    return _result("restriction", data, sequence, description,
                   parameters={"enzymes": list(enzymes), "topology": topology,
                               "include_possible": include_possible},
                   limitations=["Methylation notes are informational; blocking by methylation and star activity are not modelled."])


def design_primers(record_or_sequence: Any, target_start: int = 0,
                   target_end: Optional[int] = None, *, config: Any = None) -> AnalysisResult:
    """Primer-pair design; ``config`` is a ``ProfessionalPrimerConfig``."""
    from analyses import primer as P

    sequence, description, _ = _sequence_of(record_or_sequence)
    config = config if config is not None else P.ProfessionalPrimerConfig()
    data = P.design_primers_professional(sequence, target_start, target_end, config=config)
    thermo = config.thermodynamics
    return _result("primer", data, sequence, description,
                   parameters={"target_start": target_start, "target_end": target_end,
                               "config": config},
                   limitations=[
                       "Tm values are nearest-neighbor model predictions "
                       f"({getattr(thermo, 'nn_parameter_set', 'santalucia_1998')}, "
                       f"salt model {getattr(thermo, 'salt_model', 'monovalent_equivalent')}); "
                       "annealing temperatures must be determined experimentally.",
                   ])


def compare_sequences(reference: str, variant: str, *, reference_name: str = "Reference",
                      variant_name: str = "Variant") -> AnalysisResult:
    from analyses.mutation import CODING_FRAME_ASSUMPTION, compare_sequences as _compare

    ref = normalize_dna(reference)
    data = _compare(ref, normalize_dna(variant), reference_name=reference_name,
                    variant_name=variant_name)
    return _result("mutation", data, ref, reference_name,
                   parameters={"alignment": "global (Needleman-Wunsch), match 2, mismatch -1, gap -2"},
                   warnings=list(data.notes), limitations=[CODING_FRAME_ASSUMPTION])


def analyze_vcf(path: str, *, contig_references: Optional[Dict[str, str]] = None,
                split_multiallelic: bool = True) -> AnalysisResult:
    """Parse, normalize (when references are given) and summarize a VCF file."""
    from analyses import variant_analysis as V

    variants, metadata, samples = V.parse_vcf_file(path, split_multiallelic=split_multiallelic)
    warnings: List[str] = []
    reference_check = None
    if contig_references:
        reference_check = V.check_reference_alleles(variants, contig_references)
        if reference_check["mismatched"]:
            warnings.append(f"{reference_check['mismatched']} REF allele(s) do not match the reference.")
        variants = V.normalize_variants(variants, contig_references=contig_references)
    else:
        warnings.append("No reference sequence given: variants were not left-normalized.")
    statistics = V.calculate_variant_statistics(variants)
    data = {"variants": variants, "samples": samples, "statistics": statistics,
            "reference_check": reference_check, "header_line_count": len(metadata)}
    from core.models import sha256_file
    result = _result("variant", data, "", path,
                     parameters={"split_multiallelic": split_multiallelic,
                                 "normalized": bool(contig_references)},
                     warnings=warnings)
    result.parameters["vcf_sha256"] = sha256_file(path)
    return result


def design_crispr_guides(record_or_sequence: Any, *, system: str = "SpCas9",
                         selected_orf: Optional[Any] = None, config: Any = None) -> AnalysisResult:
    from analyses.crispr import run_crispr_analysis

    sequence, description, _ = _sequence_of(record_or_sequence)
    data = run_crispr_analysis(sequence, sequence_name=description, selected_orf=selected_orf,
                               system=system, config=config)
    return _result("crispr", data, sequence, description,
                   parameters={"system": system, "orf_restricted": selected_orf is not None},
                   limitations=["Off-target search covers only the supplied sequence, not a genome."])


__all__ = [
    "load_records", "sequence_overview", "find_orfs", "analyze_rna", "analyze_protein_from_cds",
    "restriction_digest", "design_primers", "compare_sequences", "analyze_vcf",
    "design_crispr_guides",
]
