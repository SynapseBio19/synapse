"""Terminal dialogs of analyses/primer.py (moved here in P2-2).

The scientific functions stay in ``analyses.primer``; this module draws the
screens, reads the answers and calls them.
"""

from __future__ import annotations

from analyses.primer import (  # noqa: F401
    NN_PARAMETER_SETS,
    PrimerDesignReport,
    ProfessionalPrimerConfig,
    ProfessionalPrimerDesignReport,
    RestrictionSitePrimerDesign,
    ThermodynamicConfig,
    analyze_primer,
    analyze_primer_professional,
    calculate_effective_ionic_strength,
    design_primers_for_restriction_site,
    design_primers_professional,
    normalize_dna_sequence,
    restriction_site_choices,
    summarize_primer_pair,
    summarize_professional_primer_pair,
    validate_target_region,
)
from ui import style
from ui.style import major_rule


def print_primer_report(report: PrimerDesignReport) -> None:
    """Print the complete Synapse primer report."""
    style.header('PRIMER DESIGN')
    style.field('Sequence length', f"{report.sequence_length} bp")
    style.field('Target region', f"{report.target_start + 1}-{report.target_end} (1-based; {report.target_length} bp)")
    style.section("THERMODYNAMIC SETTINGS")
    style.field('Tm model', f"{report.tm_model}")
    conditions = report.thermodynamic_conditions
    style.field('Na+ / salt', f"{conditions['salt_mM']:.2f} mM")
    style.field('Mg2+', f"{conditions['mg_mM']:.2f} mM")
    style.field('dNTPs', f"{conditions['dntp_mM']:.2f} mM")
    style.field('Primer concentration', f"{conditions['primer_concentration_nM']:.1f} nM")
    style.field('Symmetry correction', f"{'YES' if conditions['include_symmetry_correction'] else 'NO'}")
    style.section("CANDIDATE GENERATION")
    style.field('Forward evaluated', f"{report.forward_candidates_evaluated}")
    style.field('Forward accepted', f"{report.forward_candidates_accepted}")
    style.field('Reverse evaluated', f"{report.reverse_candidates_evaluated}")
    style.field('Reverse accepted', f"{report.reverse_candidates_accepted}")
    style.section("PAIR EVALUATION")
    style.field('Pairs evaluated', f"{report.pairs_evaluated}")
    style.field('Pairs accepted/ranked', f"{report.pairs_accepted}")
    print()
    if report.ranked_pairs:
        for index, pair in enumerate(report.ranked_pairs, start=1):
            print(summarize_primer_pair(pair, index))
            print()
    else:
        style.section("NO ACCEPTABLE PRIMER PAIR FOUND", blank_before=False)
        print()
        if report.failure_reasons:
            style.section("MAIN LIMITING FACTORS", blank_before=False)
            for reason in report.failure_reasons:
                print(f"- {reason}")
            print()
    if report.warnings:
        style.section("WARNINGS", blank_before=False)
        for warning in report.warnings:
            print(f"- {warning}")
    print(major_rule())


def _print_thermodynamic_conditions(
    conditions: ThermodynamicConfig,
) -> None:
    """Print the exact thermodynamic assumptions."""

    style.section('THERMODYNAMIC CONDITIONS')

    style.field('Na+', f"{conditions.sodium_mM:.3f} mM")

    style.field('Mg2+', f"{conditions.magnesium_mM:.3f} mM")

    style.field('dNTPs', f"{conditions.dntp_mM:.3f} mM")

    style.field('Primer concentration', f"{conditions.primer_nM:.3f} nM")

    effective = (
        calculate_effective_ionic_strength(
            sodium_mM=conditions.sodium_mM,
            magnesium_mM=conditions.magnesium_mM,
            dntp_mM=conditions.dntp_mM,
            salt_model=conditions.salt_model,
        )
    )

    style.field('Effective salt', f"{effective:.3f} mM")

    style.field('Concentration mode', f"{conditions.dna_concentration_mode}")

    print(
        f"Salt model:             "
        f"{conditions.salt_model}"
        + (
            " (von Ahsen et al. 2001: Na + 120*sqrt(Mg - dNTP))"
            if conditions.salt_model == "monovalent_equivalent"
            else ""
        )
    )

    style.field('NN parameters', f"{conditions.nn_parameter_set} - "
        f"{NN_PARAMETER_SETS[conditions.nn_parameter_set]['reference']}")


def _primer_ui_store_result(pipeline, result):
    """Store a primer result when the active pipeline supports it."""
    if pipeline is None:
        return

    import logging
    logger = logging.getLogger("synapse.primer")
    try:
        context = getattr(pipeline, "context", None)
        store = getattr(context, "store_primer_analysis", None)
        if callable(store):
            store(result)
        else:
            logger.warning("Pipeline has no store_primer_analysis(); primer result not stored.")
    except Exception:
        logger.exception("Could not store primer result in the pipeline.")


def _primer_ui_read_target(sequence):
    """Read an optional target region."""
    style.section('TARGET REGION')
    print("Enter 1-based, inclusive positions (the first base is 1).")
    print("Press Enter to use the complete sequence.")
    print()

    start_input = style.ask('Target start [1]')
    end_input = style.ask(f"Target end [{len(sequence)}]")

    try:
        first = 1 if not start_input else int(start_input)
        last = len(sequence) if not end_input else int(end_input)
        if not 1 <= first <= last <= len(sequence):
            raise ValueError(
                f"Positions must satisfy 1 <= start <= end <= {len(sequence)}."
            )
        # The design functions take a 0-based, half-open interval.
        target_start, target_end = first - 1, last
        validate_target_region(sequence, target_start, target_end)
        return target_start, target_end
    except (TypeError, ValueError) as exc:
        style.error(f"{exc}")
        return None


def _primer_ui_read_conditions():
    """Read thermodynamic conditions using documented defaults."""
    style.section('THERMODYNAMIC CONDITIONS')
    print("Press Enter to keep the default value.")
    print()

    try:
        sodium_input = style.ask('Na+ / salt [50.0 mM]')
        mg_input = style.ask('Mg2+ [1.5 mM]')
        dntp_input = style.ask('dNTPs [0.2 mM]')
        primer_input = style.ask('Primer concentration [250 nM]')

        sodium_mM = 50.0 if not sodium_input else float(sodium_input)
        magnesium_mM = 1.5 if not mg_input else float(mg_input)
        dntp_mM = 0.2 if not dntp_input else float(dntp_input)
        primer_nM = 250.0 if not primer_input else float(primer_input)

        conditions = ThermodynamicConfig(
            sodium_mM=sodium_mM,
            magnesium_mM=magnesium_mM,
            dntp_mM=dntp_mM,
            primer_nM=primer_nM,
        )
        conditions.validate()
        return conditions
    except (TypeError, ValueError) as exc:
        style.error(f"{exc}")
        return None


def _primer_ui_select_tm_model():
    """Select the Tm model."""
    methods = {
        "": "nearest_neighbor",  # Enter takes the [1] the prompt shows
        "1": "nearest_neighbor",
        "2": "wallace",
        "3": "salt_adjusted",
    }
    while True:
        style.menu([
            ("Tm model", [('1', 'Nearest-neighbor'), ('2', 'Wallace'), ('3', 'Salt-adjusted')]),
        ], zero="Cancel")
        choice = style.ask('Selection [1]')
        if choice == "0":
            return None
        if choice in methods:
            return methods[choice]
        style.invalid_selection()


def _primer_ui_single_analysis(pipeline=None):
    """Analyze one existing primer sequence."""
    style.header("Single primer analysis")

    primer_input = style.ask('Primer sequence')

    if not primer_input:
        style.error("Primer sequence is empty.")
        return None

    try:
        primer = normalize_dna_sequence(primer_input)
    except (TypeError, ValueError) as exc:
        style.error(f"{exc}")
        return None

    method = _primer_ui_select_tm_model()
    if method is None:
        return None

    conditions = _primer_ui_read_conditions()

    if conditions is None:
        return None

    try:
        if method == "nearest_neighbor":
            result = analyze_primer_professional(
                primer,
                conditions=conditions,
            )

            style.header('PRIMER ANALYSIS')
            style.field('Sequence', f"{result['sequence']}")
            style.field('Length', f"{result['length']} bp")
            style.field('GC content', f"{result['gc_percent']:.2f}%")
            style.field('Nearest-neighbor Tm', f"{result['nearest_neighbor_tm']:.2f} °C")
            style.field('Wallace Tm', f"{result['wallace_tm']:.2f} °C")
            style.field('Salt-adjusted Tm', f"{result['salt_adjusted_tm']:.2f} °C")
            print(f"ΔH:                     {result['delta_h_kcal_mol']:.2f} kcal/mol")
            print(f"ΔS:                     {result['delta_s_cal_mol_k']:.2f} cal/(mol·K)")
            print(f"ΔG at 25 °C:            {result['delta_g_kcal_mol']:.2f} kcal/mol")
            style.field('Effective salt', f"{result['effective_salt_mM']:.2f} mM")
            print(f"3' GC clamp:            {result['gc_clamp']} bp")
            style.field('Longest homopolymer', f"{result['max_homopolymer']}")
            style.field('Longest tandem repeat', f"{result['max_tandem_repeat_copies']}")
            style.field('Hairpin stem', f"{result['hairpin'].max_stem} bp "
                f"({result['hairpin'].severity})")
            style.field('Self-dimer stem', f"{result['self_dimer'].max_stem} bp "
                f"({result['self_dimer'].severity})")
            print(major_rule())

        else:
            result = analyze_primer(
                primer,
                tm_method=method,
                salt_mM=conditions.sodium_mM,
                mg_mM=conditions.magnesium_mM,
                dntp_mM=conditions.dntp_mM,
                primer_concentration_nM=conditions.primer_nM,
            )

            style.header('PRIMER ANALYSIS')
            style.field('Sequence', f"{result['sequence']}")
            style.field('Length', f"{result['length']} bp")
            style.field('GC content', f"{result['gc_percent']:.2f}%")
            style.field('Tm', f"{result['tm']:.2f} °C")
            style.field('Tm model', f"{result['tm_model']}")
            print(f"3' GC clamp:            {result['gc_clamp']} bp")
            style.field('Longest homopolymer', f"{result['max_homopolymer']}")
            style.field('Longest tandem repeat', f"{result['max_tandem_repeat_copies']}")
            style.field('Hairpin stem', f"{result['hairpin'].max_stem} bp "
                f"({result['hairpin'].severity})")
            style.field('Self-dimer stem', f"{result['self_dimer'].max_stem} bp "
                f"({result['self_dimer'].severity})")
            print(major_rule())

        return result

    except (TypeError, ValueError, ZeroDivisionError) as exc:
        style.error(f"{exc}")
        return None


def print_professional_primer_report(
    report: ProfessionalPrimerDesignReport,
) -> None:
    """Print the improved professional report with explicit screening status."""
    style.header('PROFESSIONAL PRIMER DESIGN')
    style.field('Sequence length', f"{report.sequence_length} bp")
    style.field('Target region', f"{report.target_start + 1}-{report.target_end} "
        f"(1-based; {report.target_length} bp)")
    style.field('Tm method', f"{report.tm_method}")
    _print_thermodynamic_conditions(report.thermodynamic_conditions)

    style.section('CANDIDATE GENERATION')
    style.field('Forward evaluated', f"{report.forward_candidates_evaluated}")
    style.field('Forward selected', f"{report.forward_candidates_accepted}")
    style.field('Reverse evaluated', f"{report.reverse_candidates_evaluated}")
    style.field('Reverse selected', f"{report.reverse_candidates_accepted}")

    style.section('SCREENING POLICY')
    print("Tier 1: strict preferred criteria")
    print("Tier 2: relaxed Tm window (+/- 2 °C)")
    print("Tier 3: best measured fallback candidates")
    print("3' GC clamp and secondary structure are scored, not blindly discarded.")

    style.section('PAIR EVALUATION')
    style.field('Pairs evaluated', f"{report.pairs_evaluated}")
    style.field('Pairs accepted/ranked', f"{report.pairs_accepted}")

    if report.ranked_pairs:
        print()
        for index, advanced_pair in enumerate(report.ranked_pairs, start=1):
            print(summarize_professional_primer_pair(advanced_pair, index))
            print()
    else:
        print()
        style.warning("No acceptable primer pair was found.")
        if report.failure_reasons:
            style.section("Main limiting factors")
            for reason in report.failure_reasons:
                print(f"- {reason}")

    if report.warnings:
        style.section("Warnings and diagnostics")
        for warning in report.warnings:
            style.warning(warning)

    print(major_rule())


def _primer_ui_design_menu_v4(dna, header, pipeline):
    """Improved compact Primer Design module menu."""
    sequence = normalize_dna_sequence(dna)

    while True:
        style.header("Primer design")
        style.field('Input', f"{header}")
        style.field('Length', f"{len(sequence)} bp")
        style.menu([
            ("Design", [('1', 'Design Primers'),
                        ('2', 'Design Primers Around a Restriction Site')]),
            ("Analysis", [('3', 'Single Primer Analysis')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            target = _primer_ui_read_target(sequence)
            if target is None:
                continue
            result = _primer_ui_design(
                sequence,
                target[0],
                target[1],
                pipeline,
            )
            if result is not None:
                return result
            continue

        if choice == "3":
            result = _primer_ui_single_analysis(pipeline)
            if result is not None:
                return result
            continue

        if choice == "2":
            result = _primer_ui_restriction_site_design(sequence, pipeline)
            if result is not None:
                return result
            continue

        style.invalid_selection()


def run_primer_analysis_ui(
    dna,
    pipeline=None,
    header="Sequence",
):
    """Module-owned Primer Design UI entry point."""
    try:
        normalize_dna_sequence(dna)
    except (TypeError, ValueError) as exc:
        style.header("Primer design")
        style.error(f"{exc}")
        return None

    return _primer_ui_design_menu_v4(
        dna,
        header,
        pipeline,
    )


def _v7_primer_ui_design(sequence, target_start, target_end, pipeline):
    """Run all three Tm models through the professional workflow."""
    tm_method = _primer_ui_select_tm_model()

    if tm_method is None:
        return None

    conditions = _primer_ui_read_conditions()
    if conditions is None:
        return None

    print()
    print("RUNNING PRIMER DESIGN...")
    print()

    try:
        config = ProfessionalPrimerConfig(
            thermodynamics=conditions,
            tm_method=tm_method,
        )

        report = design_primers_professional(
            sequence,
            target_start=target_start,
            target_end=target_end,
            config=config,
        )

        print_professional_primer_report(report)
        _primer_ui_store_result(pipeline, report)
        return report

    except (TypeError, ValueError, ZeroDivisionError) as exc:
        style.error(f"{exc}")
        return None


_primer_ui_design = _v7_primer_ui_design


def _primer_ui_restriction_site_design(sequence, pipeline):
    """Design primers around one confirmed restriction site (workflow 10.2)."""
    from data.restriction_enzymes import RESTRICTION_ENZYMES, TYPE_IIS_ENZYMES

    style.section("RESTRICTION SITE")
    known = {name.upper(): name for name in list(RESTRICTION_ENZYMES) + list(TYPE_IIS_ENZYMES)}
    answer = style.ask("Enzyme (for example EcoRI or BsaI)").strip()
    enzyme = known.get(answer.upper())
    if enzyme is None:
        style.error(f"Unknown restriction enzyme: {answer or '(none)'}.")
        return None

    choices = restriction_site_choices(sequence, enzyme)
    if not choices:
        style.error(f"{enzyme} has no confirmed site in this sequence.")
        return None

    while True:
        style.section(f"CONFIRMED {enzyme.upper()} SITES")
        style.menu([(None, [
            (str(number), f"Position {item['site_start'] + 1}-{item['site_end']}  "
                          f"{item['site_sequence']}  (cut after {item['top_cut']})")
            for number, item in enumerate(choices, 1)
        ])], zero=style.CANCEL_LABEL)
        selection = style.ask_choice()
        if selection == "0":
            return None
        if selection.isdigit() and 1 <= int(selection) <= len(choices):
            chosen = choices[int(selection) - 1]
            break
        style.invalid_selection()

    flank_input = style.ask("Bases searched on each side [300]")
    try:
        flank = 300 if not flank_input else int(flank_input)
    except ValueError:
        style.error("The flank must be a whole number of bases.")
        return None

    tm_method = _primer_ui_select_tm_model()
    if tm_method is None:
        return None
    conditions = _primer_ui_read_conditions()
    if conditions is None:
        return None

    try:
        config = ProfessionalPrimerConfig(thermodynamics=conditions, tm_method=tm_method)
        design = design_primers_for_restriction_site(
            sequence, enzyme, chosen["site_start"], flank=flank, config=config,
        )
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        style.error(f"{exc}")
        return None

    print_professional_primer_report(design.report)
    print_restriction_site_design(design)
    _primer_ui_store_result(pipeline, design)
    return design


def print_restriction_site_design(design: RestrictionSitePrimerDesign) -> None:
    """Print the site, the search regions and each amplicon's digest."""
    style.section("RESTRICTION SITE AND DIAGNOSTIC DIGEST")
    style.field("Enzyme", design.enzyme)
    style.field("Site", f"{design.site_start + 1}-{design.site_end} (1-based)")
    style.field("Forward primers in", f"{design.forward_region[0] + 1}-{design.forward_region[1]}")
    style.field("Reverse primers in", f"{design.reverse_region[0] + 1}-{design.reverse_region[1]}")
    print()
    widths = [6, 18, 10, 7]
    print(style.columns(["Rank", "Amplicon", "Length", "Sites", "Fragments (bp)"], widths))
    for item in design.amplicon_digests:
        print(style.columns([
            item["rank"],
            f"{item['amplicon_start'] + 1}-{item['amplicon_end']}",
            f"{item['amplicon_length']} bp",
            item["sites_in_amplicon"],
            " + ".join(str(length) for length in item["fragment_lengths"]),
        ], widths))
    print()
    style.notes(design.warnings, prefix=style.WARNING_PREFIX)
