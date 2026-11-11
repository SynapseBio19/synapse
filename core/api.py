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

Sequence types (audit D14): every entry point states which types it
accepts and checks its input with ``core.alphabets``. A sequence whose type
cannot be told from its letters needs ``sequence_type=``; a wrong type
raises ``core.alphabets.SequenceTypeError`` (a ``ValueError``). Nucleotide
analyses accept RNA and analyze it in the DNA alphabet (U read as T),
recording that as a warning. Restriction digests and CRISPR guide design
accept DNA only; ``analyze_protein`` takes protein sequences.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

from core.alphabets import DNA, PROTEIN, RNA, normalize_sequence, require_sequence_type
from core.models import AnalysisResult, SequenceRecord, sha256_text

_NUCLEOTIDE = (DNA, RNA)
_RNA_AS_DNA = "RNA input was analyzed in the DNA alphabet (U read as T)."


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


def _enzyme_database_version() -> str:
    from data.restriction_enzymes import ENZYME_DATABASE_VERSION
    return ENZYME_DATABASE_VERSION


def _sequence_of(record_or_sequence: Any, *, accepted: Sequence[str] = _NUCLEOTIDE,
                 sequence_type: Optional[str] = None,
                 context: str = "This analysis") -> tuple[str, str, str, List[str]]:
    """Return (sequence, description, topology, warnings) after the type check.

    A record's own declared type is used unless *sequence_type* is given.
    RNA accepted by a nucleotide analysis is converted to the DNA alphabet,
    with a warning.
    """
    if isinstance(record_or_sequence, SequenceRecord):
        raw, description, topology = (record_or_sequence.sequence, record_or_sequence.identifier,
                                      record_or_sequence.topology)
        declared = sequence_type or record_or_sequence.declared_type
    else:
        raw, description, topology = str(record_or_sequence), "sequence", "linear"
        declared = sequence_type
    result = require_sequence_type(raw, accepted, declared=declared, context=context)
    sequence = normalize_sequence(raw)
    warnings: List[str] = []
    if result.sequence_type == RNA and DNA in accepted:
        sequence = sequence.replace("U", "T")
        warnings.append(_RNA_AS_DNA)
    return sequence, description, topology, warnings


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
                      profile_step: int = 1, sequence_type: Optional[str] = None) -> AnalysisResult:
    from analyses.sequence import analyze_sequence

    sequence, description, _, warnings = _sequence_of(
        record_or_sequence, sequence_type=sequence_type, context="The sequence overview")
    data = analyze_sequence(sequence, profile_window=profile_window, profile_step=profile_step)
    return _result("sequence", data, sequence, description,
                   parameters={"profile_window": profile_window, "profile_step": profile_step},
                   warnings=warnings)


def find_orfs(record_or_sequence: Any, *, min_length: Optional[int] = None,
              start_codon_policy: Optional[str] = None, nested: Optional[bool] = None,
              ranking: Optional[str] = None, include_partial: bool = True,
              genetic_code: int = 1, top_n: int = 5,
              sequence_type: Optional[str] = None) -> AnalysisResult:
    """ORF search; ``None`` means the module default (see docs/METHODS.md)."""
    from analyses import orf as O

    sequence, description, topology, warnings = _sequence_of(
        record_or_sequence, sequence_type=sequence_type, context="ORF search")
    data = O.analyze_orfs_enhanced(
        sequence,
        include_partial=include_partial,
        min_length=O.DEFAULT_MIN_ORF_LENGTH if min_length is None else min_length,
        top_n=top_n,
        start_codon_policy=start_codon_policy,
        nested=O.DEFAULT_REPORT_NESTED if nested is None else nested,
        ranking=ranking or O.DEFAULT_RANKING,
        genetic_code=genetic_code,
        topology=topology,
    )
    limitations = ["ORFs are sequence features, not evidence of protein-coding genes."]
    if topology == "circular":
        limitations.append("Circular sequence: ORFs across the origin are included "
                           "(wraps_origin, stop_on_circle); partial ORFs are not reported.")
    return _result("orf", data, sequence, description,
                   parameters=data.get("parameters", {}), warnings=warnings,
                   limitations=limitations)


def analyze_rna(record_or_sequence: Any, *, profile_window: int = 100,
                profile_step: int = 1, sequence_type: Optional[str] = None,
                genetic_code: int = 1) -> AnalysisResult:
    """RNA analysis of a DNA template or an RNA sequence.

    ``genetic_code`` selects the NCBI table for the three-frame codon view
    (audit P1-9); pass the selected ORF's table to keep the preview
    consistent with that ORF's own translation.
    """
    from analyses.rna import analyze_rna as _analyze_rna

    sequence, description, _, warnings = _sequence_of(
        record_or_sequence, sequence_type=sequence_type, context="RNA analysis")
    data = _analyze_rna(sequence, profile_window=profile_window, profile_step=profile_step,
                        genetic_code=genetic_code)
    return _result("rna", data, sequence, description,
                   parameters={"profile_window": profile_window, "profile_step": profile_step,
                               "genetic_code": genetic_code},
                   warnings=warnings)


def analyze_protein_from_cds(cds: Any, *, initiator_start: bool = True,
                             genetic_code: Optional[int] = None,
                             description: str = "CDS",
                             sequence_type: Optional[str] = None) -> AnalysisResult:
    """Translate a coding sequence (frame +1) and analyze the protein.

    ``cds`` is a DNA or RNA string or an ORF dict from ``find_orfs``. An
    ORF is translated with its own genetic code; an explicit, different
    ``genetic_code`` raises ``ValueError``. A string uses ``genetic_code``
    (default 1). A protein sequence is refused; use ``analyze_protein``.
    """
    from analyses.orf import is_initiation_codon, orf_has_start_codon as O_has_start
    from analyses.protein import analyze_complete_protein

    if isinstance(cds, dict):
        orf_code = int(cds.get("genetic_code") or 1)
        if genetic_code is not None and int(genetic_code) != orf_code:
            raise ValueError(f"The ORF was found with NCBI table {orf_code}; "
                             f"translating it with table {genetic_code} would be inconsistent.")
        genetic_code = orf_code
        if not O_has_start(cds):
            initiator_start = False
        cds = cds["sequence"]
    genetic_code = 1 if genetic_code is None else int(genetic_code)

    sequence, _, _, warnings = _sequence_of(cds, sequence_type=sequence_type,
                                            context="Translation")
    data = analyze_complete_protein("", rna=sequence.replace("T", "U"),
                                    initiator_start=initiator_start, genetic_code=genetic_code)
    if initiator_start and len(sequence) >= 3 and not is_initiation_codon(sequence[:3], genetic_code):
        warnings.append(f"First codon {sequence[:3]} is not an initiation codon of NCBI table "
                        f"{genetic_code}; it was translated as an internal codon, not as Met.")
    return _result("protein", data, sequence, description,
                   parameters={"initiator_start": initiator_start, "genetic_code": genetic_code},
                   warnings=warnings,
                   limitations=["Physicochemical values are sequence-based estimates."])


def analyze_protein(record_or_sequence: Any, *, sequence_type: Optional[str] = None) -> AnalysisResult:
    """Analyze a protein sequence (composition, physicochemical estimates).

    A sequence of letters that are also nucleotide codes (e.g. ``MKVHDANWYSK``)
    needs ``sequence_type="PROTEIN"`` (or a record declared as protein).
    """
    from analyses.protein import analyze_complete_protein

    protein, description, _, warnings = _sequence_of(
        record_or_sequence, accepted=(PROTEIN,), sequence_type=sequence_type,
        context="Protein analysis")
    data = analyze_complete_protein(protein)
    return _result("protein", data, protein, description, warnings=warnings,
                   limitations=["Physicochemical values are sequence-based estimates."])


def restriction_digest(record_or_sequence: Any, enzymes: Sequence[str], *,
                       topology: Optional[str] = None,
                       include_possible: bool = False,
                       sequence_type: Optional[str] = None) -> AnalysisResult:
    """Restriction digest of a DNA sequence (restriction enzymes cut DNA; RNA is refused)."""
    from analyses.restriction import perform_digest

    sequence, description, record_topology, warnings = _sequence_of(
        record_or_sequence, accepted=(DNA,), sequence_type=sequence_type,
        context="A restriction digest")
    topology = topology or record_topology
    data = perform_digest(sequence, list(enzymes), topology=topology,
                          include_possible=include_possible)
    return _result("restriction", data, sequence, description,
                   parameters={"enzymes": list(enzymes), "topology": topology,
                               "include_possible": include_possible,
                               "enzyme_database_version": _enzyme_database_version()},
                   warnings=warnings,
                   limitations=["Methylation notes are informational; blocking by methylation and star activity are not modeled."])


def design_primers(record_or_sequence: Any, target_start: int = 0,
                   target_end: Optional[int] = None, *, config: Any = None,
                   sequence_type: Optional[str] = None) -> AnalysisResult:
    """Primer-pair design; ``config`` is a ``ProfessionalPrimerConfig``.

    An RNA template is designed on in the DNA alphabet (its cDNA sequence).
    """
    from analyses import primer as P

    sequence, description, _, warnings = _sequence_of(
        record_or_sequence, sequence_type=sequence_type, context="Primer design")
    config = config if config is not None else P.ProfessionalPrimerConfig()
    data = P.design_primers_professional(sequence, target_start, target_end, config=config)
    thermo = config.thermodynamics
    return _result("primer", data, sequence, description,
                   parameters={"target_start": target_start, "target_end": target_end,
                               "config": config},
                   warnings=warnings,
                   limitations=[
                       "Tm values are nearest-neighbor model predictions "
                       f"({getattr(thermo, 'nn_parameter_set', 'santalucia_1998')}, "
                       f"salt model {getattr(thermo, 'salt_model', 'monovalent_equivalent')}); "
                       "annealing temperatures must be determined experimentally.",
                   ])


def compare_sequences(reference: str, variant: str, *, reference_name: str = "Reference",
                      variant_name: str = "Variant",
                      sequence_type: Optional[str] = None) -> AnalysisResult:
    """Nucleotide comparison; *sequence_type* applies to both sequences."""
    from analyses.mutation import CODING_FRAME_ASSUMPTION, compare_sequences as _compare

    ref, _, _, ref_warnings = _sequence_of(reference, sequence_type=sequence_type,
                                           context="Sequence comparison (reference)")
    alt, _, _, alt_warnings = _sequence_of(variant, sequence_type=sequence_type,
                                           context="Sequence comparison (variant)")
    data = _compare(ref, alt, reference_name=reference_name, variant_name=variant_name)
    warnings = list(dict.fromkeys(ref_warnings + alt_warnings)) + list(data.notes)
    return _result("mutation", data, ref, reference_name,
                   parameters={"alignment": "global (Needleman-Wunsch), match 2, mismatch -1, gap -2"},
                   warnings=warnings, limitations=[CODING_FRAME_ASSUMPTION])


def analyze_vcf(path: str, *, contig_references: Optional[Dict[str, str]] = None,
                split_multiallelic: bool = True, contig_matching: str = "exact") -> AnalysisResult:
    """Parse, normalize (when references are given) and summarize a VCF file.

    ``contig_references`` maps reference contig names to sequences.
    ``contig_matching`` is "exact" (default) or "chr_prefix", which also pairs
    "1" with "chr1"; use it only for references of the same assembly
    (mitochondrial names are never paired). ``parameters["normalized"]`` is
    True only when every variant had a reference; ``data["normalization"]``
    lists the contigs that had none.
    """
    from analyses import variant_analysis as V

    V.validate_contig_matching(contig_matching)
    variants, metadata, samples = V.parse_vcf_file(path, split_multiallelic=split_multiallelic)
    warnings: List[str] = []
    reference_check = None
    normalization: Dict[str, Any] = {
        "reference_supplied": bool(contig_references),
        "contigs_not_found": [],
        "aliased_contigs": {},
        "variants_with_reference": 0,
        "variants_without_reference": len(variants),
    }
    if contig_references:
        reference_check = V.check_reference_alleles(variants, contig_references,
                                                    contig_matching=contig_matching)
        if reference_check["mismatched"]:
            warnings.append(f"{reference_check['mismatched']} REF allele(s) do not match the reference.")
        not_found = reference_check["contigs_not_found"]
        for contig in not_found:
            count = sum(v.chromosome == contig for v in variants)
            message = (f"No reference sequence for contig '{contig}': {count} variant(s) "
                       "were not REF-checked or left-aligned.")
            suggestion = reference_check["suggested_aliases"].get(contig)
            if suggestion:
                message += (f" Reference contig '{suggestion}' differs only by the 'chr' prefix; "
                            "pass contig_matching='chr_prefix' if both describe the same assembly.")
            warnings.append(message)
        for contig, reference_name in reference_check["aliased_contigs"].items():
            warnings.append(f"Contig '{contig}' was matched to reference contig '{reference_name}' "
                            "by the 'chr' prefix rule; CHROM names are unchanged.")
        variants = V.normalize_variants(variants, contig_references=contig_references,
                                        contig_matching=contig_matching)
        without_reference = sum(v.chromosome in not_found for v in variants)
        normalization.update(
            contigs_not_found=list(not_found),
            aliased_contigs=dict(reference_check["aliased_contigs"]),
            variants_with_reference=len(variants) - without_reference,
            variants_without_reference=without_reference,
        )
    else:
        warnings.append("No reference sequence given: variants were not left-normalized.")
    statistics = V.calculate_variant_statistics(variants)
    data = {"variants": variants, "samples": samples, "statistics": statistics,
            "reference_check": reference_check, "normalization": normalization,
            "header_line_count": len(metadata)}
    from core.models import sha256_file
    result = _result("variant", data, "", path,
                     parameters={"split_multiallelic": split_multiallelic,
                                 "contig_matching": contig_matching,
                                 "normalized": bool(contig_references)
                                 and normalization["variants_without_reference"] == 0},
                     warnings=warnings)
    result.parameters["vcf_sha256"] = sha256_file(path)
    return result


def design_crispr_guides(record_or_sequence: Any, *, system: str = "SpCas9",
                         selected_orf: Optional[Any] = None, config: Any = None,
                         sequence_type: Optional[str] = None) -> AnalysisResult:
    """CRISPR guide design on a DNA target (RNA is refused)."""
    from analyses.crispr import run_crispr_analysis

    sequence, description, _, warnings = _sequence_of(
        record_or_sequence, accepted=(DNA,), sequence_type=sequence_type,
        context="CRISPR guide design")
    data = run_crispr_analysis(sequence, sequence_name=description, selected_orf=selected_orf,
                               system=system, config=config)
    return _result("crispr", data, sequence, description,
                   parameters={"system": system, "orf_restricted": selected_orf is not None},
                   warnings=warnings,
                   limitations=["Off-target search covers only the supplied sequence, not a genome."])


def plan_cloning(vector: Any, insert: Any, enzymes: Optional[Sequence[str]] = None, *,
                 top_n: int = 5, vector_topology: Optional[str] = None,
                 insert_topology: Optional[str] = None,
                 sequence_type: Optional[str] = None) -> AnalysisResult:
    """Search restriction-cloning strategies for a vector and an insert (audit S14).

    Both sequences must be DNA, because restriction enzymes cut DNA. The
    full search always runs before ``top_n`` is applied, so the reported
    counts describe the whole search rather than the displayed rows.
    """
    from analyses.cloning import automatic_cloning_planner

    vector_sequence, vector_name, vector_record_topology, vector_warnings = _sequence_of(
        vector, accepted=(DNA,), sequence_type=sequence_type, context="Cloning (vector)")
    insert_sequence, insert_name, insert_record_topology, insert_warnings = _sequence_of(
        insert, accepted=(DNA,), sequence_type=sequence_type, context="Cloning (insert)")
    vector_topology = vector_topology or vector_record_topology
    insert_topology = insert_topology or insert_record_topology
    data = automatic_cloning_planner(
        vector_sequence, insert_sequence,
        enzyme_names=list(enzymes) if enzymes else None,
        top_n=top_n,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )
    return _result("cloning", data, vector_sequence, f"{vector_name} + {insert_name}",
                   parameters={"enzymes": list(enzymes) if enzymes else None, "top_n": top_n,
                               "vector_topology": vector_topology,
                               "insert_topology": insert_topology,
                               "enzyme_database_version": _enzyme_database_version()},
                   warnings=list(dict.fromkeys(vector_warnings + insert_warnings)),
                   limitations=[
                       "Strategies are sequence-level; they do not model ligation "
                       "efficiency, methylation blocking or star activity.",
                   ])


def filter_variants(variants: Sequence[Any], criteria: Any = None, **filter_fields: Any) -> AnalysisResult:
    """Apply filtering criteria to parsed variants (audit S14).

    ``criteria`` is a ``VariantFilter``; alternatively its fields may be
    passed as keywords. Variants that fail carry the reason that rejected
    them.
    """
    from analyses import variant_analysis as V

    if criteria is None:
        criteria = V.VariantFilter(**filter_fields)
    elif filter_fields:
        raise ValueError("Pass either a VariantFilter or its fields, not both.")
    data = V.filter_variants(list(variants), criteria)
    passed = len(getattr(data, "passed", []) or [])
    return _result("variant_filter", data, "", "variant filter",
                   parameters={"criteria": criteria, "input_variants": len(variants)},
                   warnings=([] if passed else ["No variant passed the filter."]))


def mutation_report(reference: str, variant: str, *, reference_name: str = "Reference",
                    variant_name: str = "Variant", include_alignment: bool = False,
                    sequence_type: Optional[str] = None) -> AnalysisResult:
    """Compare two sequences and render the mutation report text (audit S14)."""
    from analyses.mutation import build_mutation_report

    comparison = compare_sequences(reference, variant, reference_name=reference_name,
                                   variant_name=variant_name, sequence_type=sequence_type)
    report = build_mutation_report(comparison.data, include_alignment=include_alignment)
    data = {"comparison": comparison.data, "report_text": report}
    return _result("mutation_report", data, "", reference_name,
                   parameters={**comparison.parameters, "include_alignment": include_alignment},
                   warnings=list(comparison.warnings), limitations=list(comparison.limitations))


def blast_search(record_or_sequence: Any, *, program: str = "blastn", database: str = "nt",
                 engine: str = "local", sequence_name: str = "Synapse_Query",
                 sequence_type: Optional[str] = None,
                 send_to_ncbi: bool = False, **blast_options: Any) -> AnalysisResult:
    """Run a BLAST search without the interactive menu (audit S14).

    ``engine="local"`` searches an installed BLAST+ database and sends
    nothing anywhere. ``engine="remote"`` transmits the query sequence to
    NCBI (blast.ncbi.nlm.nih.gov), where it is processed on NCBI servers;
    the menu asks before doing that, so this function requires the caller
    to pass ``send_to_ncbi=True`` instead. Without it a remote search
    raises, so a script cannot submit sequences unintentionally.
    """
    from analyses.blast import run_blast_analysis

    sequence, description, _, warnings = _sequence_of(
        record_or_sequence, accepted=(DNA, RNA, PROTEIN), sequence_type=sequence_type,
        context="A BLAST search")
    engine = str(engine).lower().strip()
    if engine == "remote" and not send_to_ncbi:
        raise ValueError(
            "A remote BLAST search sends the query sequence to NCBI, where it is "
            "processed on NCBI servers. Pass send_to_ncbi=True to confirm, or use "
            "engine='local'."
        )
    if engine == "remote":
        warnings.append("The query sequence was sent to NCBI for processing.")
    run = run_blast_analysis(sequence, sequence_name=sequence_name, program=program,
                             database=database, engine=engine, **blast_options)
    return _result("blast", run, sequence, description,
                   parameters={"program": program, "database": database, "engine": engine},
                   warnings=warnings,
                   limitations=[
                       "BLAST hits show sequence similarity; they do not prove "
                       "function, orthology or identity.",
                   ])


def blast_resume(rid: str, *, contact_ncbi: bool = False, email: Optional[str] = None,
                 **client_options: Any) -> AnalysisResult:
    """Collect an NCBI BLAST search submitted earlier, by its request ID (P1-7).

    Only the request ID is sent; the query sequence is not transmitted
    again. The call still contacts NCBI, so it requires
    ``contact_ncbi=True``. NCBI keeps a request ID for a limited time, so an
    older one may have expired.
    """
    from analyses.blast import NCBIBlastClient

    if not contact_ncbi:
        raise ValueError(
            "Resuming a search contacts NCBI to fetch the stored result. Pass "
            "contact_ncbi=True to confirm. Only the request ID is sent."
        )
    client = NCBIBlastClient(email=email) if email else NCBIBlastClient()
    run = client.resume(str(rid), **client_options)
    return _result("blast", run, "", f"resumed request {rid}",
                   parameters={"rid": rid, "engine": "remote", "resumed": True},
                   warnings=["Only the request ID was sent to NCBI; the query "
                             "sequence was not transmitted again."],
                   limitations=[
                       "BLAST hits show sequence similarity; they do not prove "
                       "function, orthology or identity.",
                   ])


def literature_search(query: str, *, sources: Optional[Sequence[str]] = None,
                      max_results: int = 30, email: Optional[str] = None,
                      send_to_providers: bool = False, **query_fields: Any) -> AnalysisResult:
    """Search the configured literature providers (audit S14).

    This always makes external requests: the query text is sent to the
    selected providers. The caller must pass ``send_to_providers=True`` to
    confirm, mirroring the privacy gate the BLAST menu applies.
    """
    from analyses.literature import LiteratureQuery, run_literature_search

    if not send_to_providers:
        raise ValueError(
            "A literature search sends the query to external providers. Pass "
            "send_to_providers=True to confirm."
        )
    fields = dict(query_fields)
    if sources is not None:
        fields["sources"] = list(sources)
    literature_query = LiteratureQuery(query=query, max_results=max_results, **fields)
    data = run_literature_search(literature_query, email=email)
    return _result("literature", data, "", query,
                   parameters={"query": query, "sources": literature_query.sources,
                               "max_results": max_results},
                   warnings=["The query was sent to external literature providers."],
                   limitations=[
                       "Provider coverage, ranking and availability differ; absence "
                       "of a record is not evidence that it does not exist.",
                   ])


__all__ = [
    "load_records", "sequence_overview", "find_orfs", "analyze_rna", "analyze_protein_from_cds",
    "analyze_protein", "restriction_digest", "design_primers", "compare_sequences", "analyze_vcf",
    "design_crispr_guides", "plan_cloning", "filter_variants", "mutation_report", "blast_search",
    "blast_resume", "literature_search",
]
