"""Terminal dialogs of analyses/crispr.py (moved here in P2-2).

The scientific functions stay in ``analyses.crispr``; this module draws the
screens, reads the answers and calls them.
"""

from __future__ import annotations

from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple
import json
from analyses.crispr import (  # noqa: F401
    BaseEditingCandidate,
    CRISPRAnalysisError,
    DEFAULT_MAX_GC,
    DEFAULT_MAX_HOMOPOLYMER,
    DEFAULT_MIN_GC,
    DesignConfig,
    DesignResult,
    GuideCandidate,
    MutationTarget,
    PrimeEditingCandidate,
    _default_export_directory,
    analyze_base_editing,
    annotate_local_off_targets,
    annotate_mutation_targets,
    batch_design_guides,
    batch_summary,
    config_editing_window,
    design_guides,
    design_prime_editing_candidates,
    export_batch_summary,
    export_csv,
    export_guides_fasta,
    export_json,
    format_base_editing_report,
    format_candidate,
    format_prime_editing_report,
    format_report,
    format_system_details,
    get_base_editor,
    has_specificity_analysis,
    list_base_editors,
    list_crispr_systems,
    load_fasta,
    normalize_sequence,
    orf_target_context,
    parse_mutation_target,
    rank_guides,
    ranking_mode_label,
    run_self_test,
)
from ui import style
from ui.style import minor_rule


def print_report(result: DesignResult, *, top_n: int = 10) -> None:
    print(format_report(result, top_n=top_n))


def _input_int(prompt: str, default: Optional[int] = None, minimum: Optional[int] = None) -> int:
    while True:
        raw = input(prompt).strip()
        if not raw and default is not None:
            return default
        try:
            value = int(raw)
            if minimum is not None and value < minimum:
                raise ValueError
            return value
        except ValueError:
            suffix = f" (>= {minimum})" if minimum is not None else ""
            style.error(f"Invalid integer{suffix}.")


def _input_float(prompt: str, default: float) -> float:
    while True:
        raw = input(prompt).strip()
        if not raw:
            return default
        try:
            return float(raw)
        except ValueError:
            style.error("Invalid number.")


def _choose_system() -> str:
    systems = list_crispr_systems()
    style.menu([('Available CRISPR systems', [
        (i,
         f"{system.name:<12} PAM={system.pam:<6} Guide={system.guide_length:<2} nt | "
         f"{system.description}")
        for i, system in enumerate(systems, 1)
    ])], zero=style.CANCEL_LABEL)
    while True:
        choice = style.ask('System')
        if choice == "0":
            return ""
        try:
            index = int(choice)
            if 1 <= index <= len(systems):
                return systems[index - 1].name
        except ValueError:
            pass
        style.invalid_selection()


def _choose_editor() -> str:
    editors = list_base_editors()
    entries = []
    for i, name in enumerate(editors, 1):
        data = get_base_editor(name)
        entries.append((i, f"{name:<5} {data['conversion'][0]}>{data['conversion'][1]} window={data['window'][0]}-{data['window'][1]}"))
    style.menu([('Available computational base-editor models', entries)], zero=style.CANCEL_LABEL)
    while True:
        choice = style.ask('Editor')
        if choice == "0":
            return ""
        try:
            index = int(choice)
            if 1 <= index <= len(editors):
                return editors[index - 1]
        except ValueError:
            pass
        style.invalid_selection()


def _manual_sequence() -> Tuple[str, str]:
    name = style.ask('Sequence name') or "Manual Sequence"
    print("Paste DNA sequence. Enter an empty line to finish.")
    lines: List[str] = []
    while True:
        line = input()
        if not line.strip():
            break
        lines.append(line)
    return name, "".join(lines)


def _print_candidate_table(candidates: Sequence[GuideCandidate], limit: int = 20) -> None:
    mode = ranking_mode_label(candidates)
    style.section(f"RANKED GUIDES ({mode.upper()})")
    print("Score = transparent computational prioritization heuristic.")
    if has_specificity_analysis(candidates):
        print("Specificity component uses only explicitly supplied local reference sequences.")
    print(minor_rule())
    widths = (3, 6, 22, 5, 3, 9, 6)
    print(style.columns(("#", "ID", "Guide", "PAM", "Str", "Pos", "GC", "Score"), widths))
    print(minor_rule())
    for i, c in enumerate(candidates[:limit], 1):
        print(style.columns((
            i, c.candidate_id, c.guide, c.pam, c.strand,
            f"{c.guide_start}-{c.guide_end}", f"{c.metrics.gc_content:.1f}%",
            f"{c.rank_score:.2f}",
        ), widths))
    if len(candidates) > limit:
        print(f"Showing {limit} of {len(candidates)} candidates.")


def _choose_candidate(candidates: Sequence[GuideCandidate]) -> Optional[GuideCandidate]:
    if not candidates:
        style.note("No guide candidates available.")
        return None
    _print_candidate_table(candidates)
    raw = style.ask('Candidate number (0 cancel)')
    try:
        index = int(raw)
        if index == 0:
            return None
        if 1 <= index <= len(candidates):
            return candidates[index - 1]
    except ValueError:
        pass
    style.invalid_selection()
    return None


def _read_local_references() -> Dict[str, str]:
    print("Enter local reference sequences as name=sequence. Empty line finishes.")
    print("For PAM-aware analysis, provide the DNA context containing the guide and its PAM.")
    print("Example: target=ACCGTTGACGAACTGAAAGCTGG")
    refs: Dict[str, str] = {}
    while True:
        line = input().strip()
        if not line:
            break
        if "=" not in line:
            style.error("Expected name=sequence.")
            continue
        name, seq = line.split("=", 1)
        try:
            normalized = normalize_sequence(seq)
            if len(normalized) < 23:
                style.warning("Reference is shorter than a 20-nt SpCas9 guide plus a 3-nt PAM. "
                    "A PAM-aware hit may not be detectable.")
            refs[name.strip()] = normalized
        except CRISPRAnalysisError as exc:
            style.error(f"{exc}")
    return refs


def _load_orf_or_region(current_sequence: str) -> Tuple[str, int, int, str, str]:
    if not current_sequence:
        raise CRISPRAnalysisError("No DNA sequence selected.")
    start = _input_int("Region start (1-based): ", 1, 1)
    end = _input_int(f"Region end [{len(current_sequence)}]: ", len(current_sequence), start)
    if end > len(current_sequence):
        raise CRISPRAnalysisError("Region end exceeds sequence length.")
    return current_sequence[start - 1:end], start, end, f"Region {start}-{end}", "region"


def _cli_design() -> None:
    style.header('CRISPR GUIDE DESIGN')
    name, sequence = _manual_sequence()
    if not sequence:
        style.error("No sequence entered.")
        return
    system = _choose_system()
    if not system:
        return
    config = DesignConfig(system=system)
    config.min_gc = _input_float(f"Minimum GC [{DEFAULT_MIN_GC}]: ", DEFAULT_MIN_GC)
    config.max_gc = _input_float(f"Maximum GC [{DEFAULT_MAX_GC}]: ", DEFAULT_MAX_GC)
    config.max_homopolymer = _input_int(f"Maximum homopolymer [{DEFAULT_MAX_HOMOPOLYMER}]: ", DEFAULT_MAX_HOMOPOLYMER, 1)
    try:
        result = design_guides(sequence, sequence_name=name, system=system, config=config)
    except CRISPRAnalysisError as exc:
        style.error(f"{exc}")
        return
    print_report(result)


def crispr_menu(
    sequence: Optional[str] = None,
    sequence_name: str = "Selected Sequence",
    selected_orf: Optional[Any] = None,
) -> Optional[DesignResult]:
    """Interactive Synapse CRISPR menu. Returns the last design result (or None).

    The analysis engine remains fully featured, while the user interface groups
    closely related operations into a smaller set of professional workflows.
    Existing main.py integration is preserved.
    """
    base_sequence = normalize_sequence(sequence) if sequence else None
    current_name = sequence_name
    region_start = 1
    region_end = len(base_sequence) if base_sequence else 0
    region_name = "Selected Region"
    feature_type = "sequence"
    selected_mutations: List[MutationTarget] = []
    local_references: Dict[str, str] = {}
    last_result: Optional[DesignResult] = None
    last_base_results: List[BaseEditingCandidate] = []
    last_prime_results: List[PrimeEditingCandidate] = []
    batch_results: Dict[str, DesignResult] = {}
    batch_sequences: Dict[str, str] = {}
    config = DesignConfig()

    if selected_orf is not None:
        extracted = orf_target_context(selected_orf, base_sequence)
        if extracted is not None:
            orf_seq, orf_start, orf_name, _orf_strand = extracted
            base_sequence = orf_seq
            current_name = orf_name
            region_start = orf_start
            region_end = orf_start + len(orf_seq) - 1
            region_name = orf_name
            feature_type = "ORF"

    def sequence_context_menu() -> None:
        nonlocal base_sequence, current_name, region_start, region_end
        nonlocal region_name, feature_type, last_result

        while True:
            style.header('SEQUENCE & REGION')
            style.menu([
                (None, [('1', 'Select DNA Sequence'), ('2', 'Select ORF / Region')]),
            ])
            choice = style.ask_choice()

            if choice == "0":
                return
            if choice == "1":
                name, seq = _manual_sequence()
                try:
                    base_sequence = normalize_sequence(seq)
                except CRISPRAnalysisError as exc:
                    style.error(f"{exc}")
                    style.pause()
                    continue
                current_name = name
                region_start = 1
                region_end = len(base_sequence)
                region_name = "Selected Region"
                feature_type = "sequence"
                last_result = None
                style.field('Selected sequence', f"{current_name} ({len(base_sequence)} bp)")
            elif choice == "2":
                if not base_sequence:
                    style.note("No DNA sequence selected.")
                    style.pause()
                    continue
                try:
                    selected_seq, start, end, name, _source = _load_orf_or_region(base_sequence)
                    # start/end are 1-based within the CURRENT base sequence,
                    # which may itself be a region; keep parent coordinates.
                    offset = region_start - 1
                    base_sequence = selected_seq
                    region_start = start + offset
                    region_end = end + offset
                    region_name = f"Region {region_start}-{region_end}"
                    feature_type = "region"
                    last_result = None
                    print(f"Selected {region_name} ({len(base_sequence)} bp).")
                except CRISPRAnalysisError as exc:
                    style.error(f"{exc}")
            else:
                style.invalid_selection()

    def batch_analysis_menu() -> None:
        nonlocal batch_results, batch_sequences, current_name, base_sequence, region_start
        nonlocal region_end, region_name, feature_type, last_result

        while True:
            style.header('BATCH ANALYSIS')
            style.menu([
                (None, [
                    ('1', 'Load FASTA / Batch Analysis'),
                    ('2', 'Select Record as Current Sequence'),
                ]),
            ])
            choice = style.ask_choice()

            if choice == "0":
                return
            if choice == "1":
                path = style.ask('FASTA path')
                if not path:
                    continue
                try:
                    records = load_fasta(path)
                    batch_results = batch_design_guides(
                        records, system=config.system, config=config
                    )
                    batch_sequences = {name: seq for name, seq in records}
                    print(f"Loaded {len(records)} FASTA record(s).")
                    for row in batch_summary(batch_results):
                        print(
                            f"{row['sequence']}: {row['length_bp']} bp | "
                            f"candidates={row['candidates']} | "
                            f"high-quality={row['high_quality']}"
                        )
                except CRISPRAnalysisError as exc:
                    style.error(f"{exc}")
            elif choice == "2":
                if not batch_results:
                    style.note("No batch analysis available. Use [1] first.")
                    style.pause()
                    continue
                records = list(batch_results.keys())
                style.menu([("Available records", list(enumerate(records, start=1)))],
                           zero=style.CANCEL_LABEL)
                raw = style.ask('Record number')
                try:
                    index = int(raw)
                except ValueError:
                    style.invalid_selection()
                    continue
                if index == 0:
                    continue
                if 1 <= index <= len(records):
                    selected_name = records[index - 1]
                    result = batch_results[selected_name]
                    base_sequence = batch_sequences[selected_name]
                    current_name = selected_name
                    region_start = result.region_start
                    region_end = result.region_end
                    region_name = result.region_name
                    feature_type = result.feature_type
                    last_result = result
                    style.field('Current sequence', f"{current_name}")
                else:
                    style.invalid_selection()
            else:
                style.invalid_selection()

    def system_menu() -> None:
        nonlocal config
        while True:
            style.header('CRISPR SYSTEM')
            style.field('Current system', f"{config.system}")
            style.menu([
                (None, [('1', 'Select CRISPR System'), ('2', 'Show System Details')]),
            ])
            choice = style.ask_choice()

            if choice == "0":
                return
            if choice == "1":
                system = _choose_system()
                if system:
                    config.system = system
                    style.field('Selected system', f"{system}")
            elif choice == "2":
                print()
                print(format_system_details(config.system))
            else:
                style.invalid_selection()

    def guide_analysis_menu() -> None:
        nonlocal last_result
        if not last_result:
            style.note("Run [4] Target Discovery first.")
            return
        while True:
            style.header('GUIDE ANALYSIS & RANKING')
            style.menu([
                (None, [('1', 'Analyze Individual Guide'), ('2', 'Rank Guides')]),
            ])
            choice = style.ask_choice()
            if choice == "0":
                return
            if choice == "1":
                candidate = _choose_candidate(last_result.candidates)
                if candidate:
                    print("\n" + format_candidate(candidate, 1))
            elif choice == "2":
                last_result.candidates = rank_guides(last_result.candidates, config=config)
                _print_candidate_table(last_result.candidates)
            else:
                style.invalid_selection()

    def report_export_menu() -> None:
        if not last_result and not batch_results:
            style.note("No analysis results available.")
            return
        while True:
            style.header('REPORTS & EXPORT')
            style.menu([
                (None, [
                    ('1', 'Full Guide Report'),
                    ('2', 'Export Guides'),
                    ('3', 'Export Batch Summary'),
                ]),
            ])
            choice = style.ask_choice()

            if choice == "0":
                return
            if choice == "1":
                if not last_result:
                    style.note("Run [4] Target Discovery first.")
                    style.pause()
                    continue
                print_report(last_result, top_n=10)
                if last_base_results:
                    print(format_base_editing_report(last_base_results))
                if last_prime_results:
                    print(format_prime_editing_report(last_prime_results, top_n=10))
            elif choice == "2":
                if not last_result:
                    style.note("Run [4] Target Discovery first.")
                    style.pause()
                    continue
                default_dir = _default_export_directory()
                style.note(f"Relative paths are saved to: {default_dir}")
                path = style.ask('Export path (.json, .csv or .fasta)')
                if not path:
                    continue
                try:
                    lower = path.lower()
                    if lower.endswith(".json"):
                        exported = export_json(last_result, path)
                    elif lower.endswith(".csv"):
                        exported = export_csv(last_result, path)
                    elif lower.endswith(".fasta") or lower.endswith(".fa"):
                        exported = export_guides_fasta(last_result, path)
                    else:
                        style.error("Use .json, .csv or .fasta/.fa extension.")
                        continue
                    style.success(f"Exported: {exported}")
                except PermissionError:
                    style.error("Export failed: permission denied.")
                    print("Choose a writable Desktop/Documents location or use a relative filename.")
                except OSError as exc:
                    style.error(f"Export failed: {exc}")
            elif choice == "3":
                if not batch_results:
                    style.note("No batch analysis available. Use [2] -> [1] Load FASTA / Batch Analysis first.")
                    style.pause()
                    continue
                default_dir = _default_export_directory()
                style.note(f"Relative paths are saved to: {default_dir}")
                path = style.ask('Batch CSV export path')
                if not path:
                    continue
                try:
                    exported = export_batch_summary(batch_results, path)
                    style.success(f"Exported: {exported}")
                except PermissionError:
                    style.error("Export failed: permission denied.")
                    print("Choose a writable Desktop/Documents location or use a relative filename.")
                except OSError as exc:
                    style.error(f"Export failed: {exc}")
            else:
                style.invalid_selection()

    while True:
        style.header('CRISPR GUIDE DESIGN')
        style.field('Sequence', f"{current_name}")
        style.field('Length', f"{len(base_sequence) if base_sequence else 0} bp")
        style.field('Region', f"{region_name} ({region_start}-{region_end})" if base_sequence else "none")
        style.field('System', f"{config.system}")
        style.menu([
            ('Input / Context', [
                ('1', 'Sequence & Region'),
                ('2', 'Batch Analysis'),
                ('3', 'CRISPR System'),
            ]),
            ('Guide Design', [('4', 'Target Discovery'), ('5', 'Guide Analysis & Ranking')]),
            ('Editing', [
                ('6', 'Mutation / Variant Targeting'),
                ('7', 'Base Editing'),
                ('8', 'Prime Editing'),
            ]),
            ('Specificity', [('9', 'Off-Target Analysis'), ('10', 'Specificity Profile')]),
            ('Output', [('11', 'Reports & Export')]),
            ('System', [('12', 'CRISPR Self-Test'), ('13', 'Supported Systems')]),
        ])
        choice = style.ask_choice()

        if choice == "0":
            # Returned to main.py so the result can be stored in the pipeline.
            return last_result

        if choice == "1":
            sequence_context_menu()

        elif choice == "2":
            batch_analysis_menu()

        elif choice == "3":
            system_menu()

        elif choice == "4":
            if not base_sequence:
                style.note("No DNA sequence selected.")
                style.pause()
                continue
            config.min_gc = _input_float(f"Minimum GC [{config.min_gc}]: ", config.min_gc)
            config.max_gc = _input_float(f"Maximum GC [{config.max_gc}]: ", config.max_gc)
            try:
                last_result = design_guides(
                    base_sequence,
                    sequence_name=current_name,
                    system=config.system,
                    region_start=region_start,
                    region_end=region_end,
                    region_name=region_name,
                    feature_type=feature_type,
                    feature_name=region_name,
                    config=config,
                    mutations=selected_mutations,
                    local_reference_sequences=local_references or None,
                )
                print(f"Found {last_result.candidate_count} candidate(s).")
                _print_candidate_table(last_result.candidates)
            except CRISPRAnalysisError as exc:
                style.error(f"{exc}")

        elif choice == "5":
            guide_analysis_menu()

        elif choice == "6":
            raw = style.ask('Variant (chr:position REF>ALT; blank cancels)')
            if not raw:
                continue
            try:
                mutation = parse_mutation_target(raw)
                selected_mutations.append(mutation)
                style.field('Added variant', f"{mutation.chromosome}:{mutation.position} "
                    f"{mutation.reference}>{mutation.alternative}")
                if last_result:
                    annotate_mutation_targets(
                        last_result.candidates,
                        selected_mutations,
                        window=config.mutation_window,
                    )
                    last_result.candidates = rank_guides(last_result.candidates, config=config)
                    _print_candidate_table(last_result.candidates)
            except CRISPRAnalysisError as exc:
                style.error(f"{exc}")

        elif choice == "7":
            if not last_result:
                style.note("Run [4] Target Discovery first.")
                style.pause()
                continue
            if not selected_mutations:
                style.note("No mutation/variant selected. Use [6] first.")
                style.pause()
                continue
            editor = _choose_editor()
            if not editor:
                continue
            last_base_results = []
            for mutation in selected_mutations:
                last_base_results.extend(
                    analyze_base_editing(last_result.candidates, mutation, editor=editor,
                                         window=config_editing_window(config))
                )
            print(format_base_editing_report(last_base_results))

        elif choice == "8":
            if not last_result or not base_sequence:
                style.note("Run [4] Target Discovery first.")
                style.pause()
                continue
            if not selected_mutations:
                style.note("No mutation/variant selected. Use [6] first.")
                style.pause()
                continue
            last_prime_results = []
            for mutation in selected_mutations:
                last_prime_results.extend(
                    design_prime_editing_candidates(
                        base_sequence,
                        mutation,
                        last_result.candidates,
                        pbs_lengths=config.prime_pbs_lengths,
                        rtt_lengths=config.prime_rtt_lengths,
                        sequence_start=region_start,
                    )
                )
            print(format_prime_editing_report(last_prime_results))

        elif choice == "9":
            if not last_result:
                style.note("Run [4] Target Discovery first.")
                style.pause()
                continue
            local_references = _read_local_references()
            if local_references:
                annotate_local_off_targets(
                    last_result.candidates,
                    local_references,
                    max_mismatches=config.max_local_off_target_mismatches,
                    max_hits=config.max_local_off_target_hits,
                    seed_length=config.seed_length,
                    require_pam=config.require_pam_for_off_targets,
                )
                last_result.candidates = rank_guides(last_result.candidates, config=config)
                total_hits = sum(len(c.local_off_target_hits) for c in last_result.candidates)
                style.success("PAM-aware local specificity analysis completed.")
                if total_hits:
                    print(f"Local PAM-compatible hits found: {total_hits}")
                    print("Top local matches:")
                    shown = 0
                    for candidate in last_result.candidates:
                        for hit in candidate.local_off_target_hits:
                            print(
                                f"  {candidate.candidate_id}: {hit['sequence']} "
                                f"strand={hit['strand']} guide={hit['guide_start']}-{hit['guide_end']} "
                                f"PAM={hit['pam']} mismatches={hit['mismatches']} "
                                f"seed_mismatches={hit['seed_mismatches']} "
                                f"classification={hit['classification']}"
                            )
                            shown += 1
                            if shown >= 20:
                                break
                        if shown >= 20:
                            break
                    if total_hits > shown:
                        print(f"Showing {shown} of {total_hits} local hits.")
                else:
                    print("No PAM-compatible local hits were found in the supplied references.")
                    print(
                        "For a positive control, provide reference DNA containing the guide plus its PAM "
                        "(for example: target=ACCGTTGACGAACTGAAAGCTGG for SpCas9/G0004)."
                    )
                print("Scientific limitation: this is not a genome-wide off-target assessment.")

        elif choice == "10":
            if not last_result:
                style.note("Run [4] Target Discovery first.")
                style.pause()
                continue
            style.section('SPECIFICITY PROFILE')
            if has_specificity_analysis(last_result.candidates):
                print("Candidates with explicit local specificity analysis are marked in the ranking.")
                print("Select the candidate to inspect its local specificity profile.")
            else:
                style.note("No local specificity analysis has been performed yet.")
                style.note("Use [9] Off-Target Analysis first and supply reference sequences.")
            candidate = _choose_candidate(last_result.candidates)
            if candidate:
                profile = candidate.specificity_profile or {
                    "analyzed": False,
                    "scope": "local supplied sequences",
                    "query": candidate.candidate_id,
                    "total_hits": 0,
                    "genome_wide": False,
                    "limitation": "No explicit reference sequence set has been searched.",
                }
                print(json.dumps(profile, indent=2))

        elif choice == "11":
            report_export_menu()

        elif choice == "12":
            run_self_test(verbose=True)

        elif choice == "13":
            style.section('SUPPORTED CRISPR SYSTEMS')
            for system in list_crispr_systems():
                print(
                    f"{system.name}: family={system.enzyme_family}, "
                    f"PAM={system.pam}, guide={system.guide_length} nt, "
                    f"PAM side={system.pam_side}"
                )

        else:
            style.invalid_selection()
