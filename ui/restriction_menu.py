"""Terminal dialogs of analyses/restriction.py (moved here in P2-2).

The scientific functions stay in ``analyses.restriction``; this module draws the
screens, reads the answers and calls them.
"""

from data.restriction_enzymes import RESTRICTION_ENZYMES
from data.restriction_enzymes import TYPE_IIS_ENZYMES
from data.restriction_enzymes import get_restriction_enzyme
from analyses.restriction import (  # noqa: F401
    analyze_restriction,
    get_digest_prediction,
    get_gel_band_prediction,
    get_restriction_map,
)
from ui import style
from ui.style import major_rule, minor_rule


def print_restriction_report(analysis):
    """Print a human-readable report without changing analysis data."""
    if not analysis:
        style.note("No restriction analysis available.")
        return

    style.header('RESTRICTION ANALYSIS')

    style.section('Sequence')
    style.field('Length', f"{analysis['sequence_length']} bp")
    style.field('Topology', f"{analysis['topology'].capitalize()}")
    style.field('Enzymes analyzed', f"{len(analysis['enzymes'])}")

    status = analysis.get("sequence_status", {})
    ambiguous = status.get("ambiguous_positions", [])
    if ambiguous:
        style.field('IUPAC ambiguous bases', f"{len(ambiguous)} "
            f"(positions: {', '.join(map(str, ambiguous[:20]))}"
            f"{'…' if len(ambiguous) > 20 else ''})")

    style.section('Restriction Enzymes')
    widths = (10, 13, 6, 10, 9)
    print(style.columns(("Enzyme", "Recognition", "Sites", "Confirmed", "Possible", "Type"), widths))
    print(minor_rule())

    for enzyme_name in analysis["enzymes"]:
        result = analysis["enzyme_results"][enzyme_name]
        print(style.columns((
            enzyme_name,
            result["recognition_sequence"],
            result["site_count"],
            result["confirmed_site_count"],
            result["possible_site_count"],
            result["cut_type"],
        ), widths))

    style.section('Restriction Map')
    sites = analysis.get("restriction_map", [])

    if not sites:
        print("No restriction sites detected.")
    else:
        for site in sites:
            overhang = site.get("overhang_sequence") or "-"
            status_text = site.get("match_status", "UNKNOWN")
            print(
                f"{site['enzyme']:<12} "
                f"site {site['site_position']:>7} "
                f"| cut {site['top_cut']:>7}/{site['bottom_cut']:<7} "
                f"| {site['cut_type']:<6} "
                f"| {site['overhang_type']:<12} "
                f"| {overhang:<8} "
                f"| {status_text}"
            )

    digest = analysis["digest"]

    style.section('Digest')
    style.field('Enzyme cut events', f"{digest['cut_count']}")
    style.field('Physical cut boundaries', f"{digest['unique_cut_count']}")
    style.field('Expected fragments', f"{digest['fragment_count']}")

    if digest["fragments"]:
        for fragment in digest["fragments"]:
            enzymes = ", ".join(fragment.get("enzymes", [])) or "-"
            print(
                f"Fragment {fragment['fragment']:>3}: "
                f"{fragment['length']:>8} bp"
                f" | Cut by: {enzymes}"
            )

    style.section('Predicted Gel Bands')
    bands = analysis.get("gel_bands", [])
    if not bands:
        style.note("No bands available.")
    else:
        for band in bands:
            print(
                f"{band['size']:>8} bp | "
                f"{band['copies']} fragment"
                f"{'s' if band['copies'] != 1 else ''}"
            )

    style.section('Unique Cutters')
    print(", ".join(analysis["unique_cutters"]) or "None detected.")

    style.section('Non-Cutting Enzymes')
    print(", ".join(analysis["non_cutters"]) or "None detected.")

    style.section('Multiple-Site Warnings')
    warnings = analysis.get("internal_sites", [])
    if not warnings:
        print("No multiple-site warnings detected.")
    else:
        for warning in warnings:
            positions = ", ".join(map(str, warning["positions"]))
            style.warning(f"{warning['enzyme']} has "
                f"{warning['site_count']} sites at {positions}")

    style.section('Recommended Enzyme Pairs')
    recommendations = analysis.get("recommended_pairs", [])
    if not recommendations:
        print("No suitable enzyme pairs detected.")
    else:
        for index, recommendation in enumerate(recommendations, 1):
            print(
                f"{index}. {recommendation['enzyme_a']} + "
                f"{recommendation['enzyme_b']} "
                f"({recommendation['score']}/{recommendation['maximum']})"
            )
            for reason in recommendation["reasons"]:
                print(f"   + {reason}")
            for warning in recommendation["warnings"]:
                print(f"   ! {warning}")

    style.section('Interpretation')
    if digest["cut_count"] == 0:
        print("No selected restriction enzyme cuts the sequence.")
    elif analysis["topology"] == "circular":
        print(
            "Circular DNA: the number of fragments equals the number "
            "of distinct physical cut boundaries."
        )
    else:
        print(
            "Linear DNA: n distinct internal cut boundaries produce "
            "n + 1 fragments."
        )

    if analysis.get("site_statistics"):
        stats = analysis["site_statistics"]
        style.field('Site events', f"{stats['total_site_events']} | "
            f"confirmed: {stats['confirmed_site_events']} | "
            f"possible: {stats['possible_site_events']}")

    print(
        "\nNote: This is sequence-level in-silico analysis. "
        "Experimental digestion depends on enzyme-specific conditions, "
        "DNA context, methylation and other laboratory factors."
    )
    print(major_rule())


def _synapse_ui_pause():
    style.pause()


def _synapse_ui_select_enzymes(
    allow_all=True,
    minimum=1,
    maximum=None,
    include_type_iis=False,
):
    enzyme_names = list(RESTRICTION_ENZYMES.keys())
    type_iis_names = list(TYPE_IIS_ENZYMES) if include_type_iis else []

    style.section("Available restriction enzymes")

    for number, enzyme_name in enumerate(enzyme_names, start=1):
        enzyme = RESTRICTION_ENZYMES[enzyme_name]
        print(
            f"{number:>2}. {enzyme_name:<10} "
            f"{enzyme['recognition_sequence']}"
        )

    if type_iis_names:
        print("\nType IIS (cut outside the site; both strands are searched):")
        for number, enzyme_name in enumerate(type_iis_names, start=len(enzyme_names) + 1):
            enzyme = TYPE_IIS_ENZYMES[enzyme_name]
            size = len(enzyme["recognition_sequence"])
            print(
                f"{number:>2}. {enzyme_name:<10} "
                f"{enzyme['recognition_sequence']}"
                f"({enzyme['cut_top'] - size}/{enzyme['cut_bottom'] - size})"
            )

    print("\nEnter enzyme names separated by commas.")
    print("Example: EcoRI,HindIII")
    if allow_all:
        print("Enter ALL to analyze all enzymes.")
    print("Enter 0 to go back.")

    while True:
        print()
        raw = style.ask("Enzymes")

        if raw == "0":
            return None

        if allow_all and raw.upper() == "ALL":
            return enzyme_names

        selected = list(dict.fromkeys(
            name.strip()
            for name in raw.split(",")
            if name.strip()
        ))

        invalid = [
            name for name in selected
            if name not in RESTRICTION_ENZYMES
            and name not in type_iis_names
        ]

        if invalid:
            style.error("Unknown enzyme(s): " + ", ".join(invalid))
            continue

        if len(selected) < minimum:
            style.error(f"Please select at least {minimum} enzyme(s).")
            continue

        if maximum is not None and len(selected) > maximum:
            style.error(f"Please select no more than {maximum} enzyme(s).")
            continue

        return selected


def _synapse_ui_select_topology(default="linear"):
    style.menu([("DNA topology", [
        (1, "Linear"),
        (2, "Circular"),
        ("Enter", f"Keep current ({str(default).capitalize()})"),
    ])])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None
        if choice == "":
            return str(default).lower()
        if choice == "1":
            return "linear"
        if choice == "2":
            return "circular"

        style.invalid_selection()


def _restriction_overview_ui(dna, pipeline=None, topology="linear"):
    style.header('RESTRICTION OVERVIEW')

    enzymes = _synapse_ui_select_enzymes()
    if enzymes is None:
        return None

    selected_topology = _synapse_ui_select_topology(topology)
    if selected_topology is None:
        return None

    try:
        analysis = analyze_restriction(
            dna,
            enzymes,
            selected_topology,
        )

        if pipeline is not None:
            try:
                pipeline.context.store_restriction_analysis(analysis)
            except AttributeError:
                pass

        print_restriction_report(analysis)
        _synapse_ui_pause()
        return analysis

    except (ValueError, TypeError, KeyError) as error:
        print()
        style.error(f"Restriction analysis failed: {error}")
        _synapse_ui_pause()
        return None


def _restriction_map_ui(dna, topology="linear"):
    style.header('RESTRICTION MAP')

    enzymes = _synapse_ui_select_enzymes()
    if enzymes is None:
        return None

    selected_topology = _synapse_ui_select_topology(topology)
    if selected_topology is None:
        return None

    try:
        sites = get_restriction_map(
            dna,
            enzyme_names=enzymes,
            topology=selected_topology,
        )

        print()
        print(minor_rule())

        if not sites:
            print("No restriction sites detected.")
        else:
            for site in sites:
                position = site.get(
                    "site_position",
                    site.get("position", "-"),
                )
                print(
                    f"{site.get('enzyme', '-'):<10} "
                    f"site {position:>6} | "
                    f"cut {site.get('top_cut', '-')}/"
                    f"{site.get('bottom_cut', '-')} | "
                    f"{site.get('cut_type', '-'):<6} | "
                    f"{site.get('overhang_sequence') or '-'}"
                )

        _synapse_ui_pause()
        return sites

    except (ValueError, TypeError, KeyError) as error:
        print()
        style.error(f"Restriction map failed: {error}")
        _synapse_ui_pause()
        return None


def _restriction_digest_gel_ui(dna, topology="linear"):
    style.header('DIGEST & GEL ANALYSIS')

    enzymes = _synapse_ui_select_enzymes(
        allow_all=False,
        minimum=1,
        include_type_iis=True,
    )
    if enzymes is None:
        return None

    selected_topology = _synapse_ui_select_topology(topology)
    if selected_topology is None:
        return None

    try:
        digest = get_digest_prediction(
            dna,
            enzymes,
            topology=selected_topology,
        )
        bands = get_gel_band_prediction(digest)

        print()
        style.field('Topology', f"{selected_topology.capitalize()}")
        style.field('Enzymes', f"{', '.join(enzymes)}")
        style.field('Cut events', f"{digest['cut_count']}")
        style.field('Cut boundaries', f"{digest['unique_cut_count']}")
        style.field('Fragments', f"{digest['fragment_count']}")

        possible = digest.get("possible_cuts", [])
        if possible:
            print()
            style.note(f"{len(possible)} ambiguity-dependent site(s) "
                       "(IUPAC symbols such as N inside the recognition "
                       "sequence) were NOT counted as cuts.")

        for enzyme_name in enzymes:
            enzyme = get_restriction_enzyme(enzyme_name) or {}
            note = enzyme.get("methylation_note")
            if note:
                print(f"Methylation note ({enzyme_name}): {note}")

        type_iis_cuts = [cut for cut in digest["cuts"] if "strand" in cut]
        if type_iis_cuts:
            style.section('Type IIS cuts')
            widths = (8, 7, 10, 10, 13)
            print(style.columns(("Enzyme", "Strand", "Site", "Top cut", "Bottom cut", "Overhang"), widths))
            for cut in type_iis_cuts:
                print(style.columns((
                    cut["enzyme"], cut["strand"], cut["site_position_1based"],
                    f"after {cut['top_cut_0based']}", f"after {cut['bottom_cut_0based']}",
                    cut["cut_info"]["overhang"] or "-",
                ), widths))
            style.note("Positions are 1-based; a cut lies after the base given. The overhang "
                       "is the top-strand sequence between the two cuts.")

        style.section('Fragments')
        for fragment in digest.get("fragments", []):
            print(
                f"Fragment {fragment['fragment']:>2}: "
                f"{fragment['length']:>8} bp"
            )

        style.section('Predicted Gel Bands')
        if not bands:
            style.note("No bands available.")
        else:
            for band in bands:
                copies = band["copies"]
                print(
                    f"{band['size']:>8} bp | "
                    f"{copies} fragment"
                    f"{'s' if copies != 1 else ''}"
                )

        _synapse_ui_pause()
        return {"digest": digest, "gel_bands": bands}

    except (ValueError, TypeError, KeyError) as error:
        print()
        style.error(f"Digest analysis failed: {error}")
        _synapse_ui_pause()
        return None


def _restriction_complete_report_ui(
    dna,
    pipeline=None,
    topology="linear",
):
    style.header('COMPLETE RESTRICTION REPORT')

    enzymes = _synapse_ui_select_enzymes()
    if enzymes is None:
        return None

    selected_topology = _synapse_ui_select_topology(topology)
    if selected_topology is None:
        return None

    try:
        analysis = analyze_restriction(
            dna,
            enzymes,
            selected_topology,
        )

        if pipeline is not None:
            try:
                pipeline.context.store_restriction_analysis(analysis)
            except AttributeError:
                pass

        print_restriction_report(analysis)
        _synapse_ui_pause()
        return analysis

    except (ValueError, TypeError, KeyError) as error:
        print()
        style.error(f"Restriction analysis failed: {error}")
        _synapse_ui_pause()
        return None


def run_restriction_analysis_ui(
    dna,
    pipeline=None,
    header="Sequence",
    topology="linear",
):
    """Small user-facing Restriction menu."""

    while True:
        style.header('RESTRICTION ANALYSIS')
        style.field("Sequence", header)
        style.field("Length", f"{len(dna)} bp")
        style.field("Topology", topology)

        style.menu([
            ("Analysis", [
                ('1', 'Restriction Overview'),
                ('2', 'Restriction Map'),
                ('3', 'Digest & Gel Analysis'),
            ]),
            ("Report", [('4', 'Complete Restriction Report')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            _restriction_overview_ui(
                dna,
                pipeline=pipeline,
                topology=topology,
            )
        elif choice == "2":
            _restriction_map_ui(dna, topology=topology)
        elif choice == "3":
            _restriction_digest_gel_ui(
                dna,
                topology=topology,
            )
        elif choice == "4":
            _restriction_complete_report_ui(
                dna,
                pipeline=pipeline,
                topology=topology,
            )
        else:
            print()
            style.invalid_selection()


def run_restriction_cloning_ui(
    dna,
    header,
    sequences,
    pipeline=None,
    topology="linear",
):
    """
    Top-level module 5 UI.

    First choose Restriction or Cloning.
    Cloning remains owned by the cloning module.
    """

    try:
        from ui.cloning_menu import run_cloning_analysis_ui
    except ImportError:
        run_cloning_analysis_ui = None

    while True:
        style.header('RESTRICTION & CLONING ANALYSIS')
        style.field("Sequence", header)
        style.field("Length", f"{len(dna)} bp")
        style.field("Topology", topology)

        style.menu([
            ("Analysis", [('1', 'Restriction Analysis'), ('2', 'Cloning Analysis')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            run_restriction_analysis_ui(
                dna=dna,
                pipeline=pipeline,
                header=header,
                topology=topology,
            )
            continue

        if choice == "2":
            if run_cloning_analysis_ui is None:
                print()
                style.error("Cloning module could not be imported.")
                _synapse_ui_pause()
                continue

            run_cloning_analysis_ui(
                vector_sequence=dna,
                vector_header=header,
                sequences=sequences,
                pipeline=pipeline,
                vector_topology=topology,
            )
            continue

        print()
        style.invalid_selection()
