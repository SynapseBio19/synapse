"""Terminal dialogs of analyses/cloning.py (moved here in P2-2).

The scientific functions stay in ``analyses.cloning``; this module draws the
screens, reads the answers and calls them.
"""

from data.restriction_enzymes import RESTRICTION_ENZYMES
from analyses.cloning import (  # noqa: F401
    _normalize_topology,
    analyze_cloning,
    automatic_cloning_planner,
    evaluate_cloning_strategy,
    get_design_reasons,
    normalize_sequence,
    validate_dna_sequence,
)
from ui import style
from ui.style import major_rule, minor_rule


def _yes_no(value):
    return "YES" if value else "NO"


def print_design_reasons(
    analysis,
):
    evaluation = get_design_reasons(
        analysis
    )

    style.section('DESIGN EVALUATION')

    for reason in evaluation[
        "positive"
    ]:
        print(
            f"+ {reason}"
        )

    for warning in evaluation[
        "warnings"
    ]:
        print(
            f"! {warning}"
        )


def ascii_orientation_visualization(
    analysis,
):
    if not analysis:
        return

    vector = analysis[
        "vector"
    ]
    insert = analysis[
        "insert"
    ]
    strategy = analysis[
        "strategy"
    ]

    style.section('CLONING ARCHITECTURE')
    print(
        "VECTOR:"
    )
    print(
        f"{vector['left']['enzyme']}"
        f" ========================== "
        f"{vector['right']['enzyme']}"
    )
    print(
        "                 | INSERT |"
    )
    print(
        "INSERT:"
    )
    print(
        f"{insert['left']['enzyme']}"
        f" >>>>>>>>>>>>>>>>>>>>>>>>>>> "
        f"{insert['right']['enzyme']}"
    )
    print()
    style.field('Orientation', f"{strategy['orientation']}")


def print_cloning_report(
    analysis,
):
    style.header('CLONING ANALYSIS')

    if not analysis:
        print()
        style.note("No analysis available.")
        print(
            major_rule()
        )
        return

    if analysis.get(
        "planner"
    ) is not None:
        planner = analysis[
            "planner"
        ]
        best = planner.get(
            "best"
        )

        if best is not None:
            _print_single_strategy_report(
                best
            )
        else:
            print_planner_report(
                planner
            )

        return

    if "vector" not in analysis:
        print()
        style.field('Status', f"{analysis.get('status', 'UNKNOWN')}")
        print(
            major_rule()
        )
        return

    if analysis.get("mode") == "single-enzyme":
        _print_single_enzyme_report(
            analysis
        )
        return

    _print_single_strategy_report(
        analysis
    )


def _pass_fail(value):
    return "PASS" if value else "FAIL"


def _print_single_enzyme_report(
    analysis,
):
    """Report for analyze_cloning(): one enzyme opens the vector and cuts the insert.

    The two-enzyme report reads left and right ends, which this result does
    not have; it shows the fields analyze_cloning() returns.
    """
    validity = analysis["validity"]
    strategy = analysis["strategy"]
    compatibility = analysis["compatibility"]
    score = analysis["score"]

    for title, site, molecule in (
        ("VECTOR", analysis["vector_site"], analysis["vector"]),
        ("INSERT", analysis["insert_site"], analysis["insert"]),
    ):
        style.section(title)
        style.field('Enzyme', f"{site['enzyme']}")
        style.field('Sites found', f"{molecule['site_count']}")
        style.field('First site position', f"{site['position'] if site['position'] is not None else '-'}")
        style.field('Recognition', f"{site['recognition_sequence']}")
        style.field('Overhang', f"{site['overhang'] or '-'} ({site['overhang_type']})")

    style.section('CLONING JUNCTION')
    style.field('Status', f"{compatibility['status']}")
    style.field('Ends', f"{compatibility['vector_end'] or '-'} -> {compatibility['insert_end'] or '-'}")
    style.field('Mechanism', f"{compatibility['mechanism']}")

    style.section('ORIENTATION')
    style.field('Strategy', f"{strategy['strategy']}")
    style.field('Orientation', f"{strategy['orientation']}")
    style.field('Directional', f"{_yes_no(strategy['directional'])}")
    style.field('Self-ligation risk', f"{'POSSIBLE' if strategy['self_ligation_risk'] else 'REDUCED'}")

    style.section('DESIGN SCORE')
    style.field('Score', f"{score['score']}/{score['maximum']}")
    style.field('Category', f"{score['category']}")
    for reason in score["reasons"]:
        print(f"- {reason}")

    style.section('VALIDATION')
    style.field('Required sites', _pass_fail(validity['all_required_sites']))
    style.field('Vector site unique', _pass_fail(validity['vector_unique']))
    style.field('Insert site unique', _pass_fail(validity['insert_unique']))
    style.field('Junction', _pass_fail(validity['junctions_compatible']))
    style.field('Overall valid', f"{'YES' if validity['valid'] else 'NO'}")

    if analysis.get("warnings"):
        style.section('WARNINGS')
        for warning in analysis["warnings"]:
            print(f"- {warning}")

    style.section('INTERPRETATION')
    if validity["valid"]:
        print("The selected design passes the implemented sequence-level cloning checks.")
    else:
        print("The selected design does not pass all implemented sequence-level cloning checks.")
    print(
        "\nThis is a sequence-level design assessment. "
        "It does not guarantee experimental cleavage, ligation, "
        "enzyme activity, methylation behavior, buffer compatibility, "
        "or wet-lab success."
    )
    print(major_rule())


def _print_single_strategy_report(
    analysis,
):
    vector = analysis[
        "vector"
    ]
    insert = analysis[
        "insert"
    ]
    strategy = analysis[
        "strategy"
    ]
    score = analysis[
        "score"
    ]
    validity = analysis[
        "validity"
    ]

    style.section('VECTOR')
    style.field('Left enzyme', f"{vector['left']['enzyme']}")
    style.field('Left position', f"{vector['left']['position']}")
    style.field('Left recognition', f"{vector['left']['recognition_sequence']}")
    style.field('Left overhang', f"{vector['left']['overhang'] or '-'}")
    style.field('Right enzyme', f"{vector['right']['enzyme']}")
    style.field('Right position', f"{vector['right']['position']}")
    style.field('Right recognition', f"{vector['right']['recognition_sequence']}")
    style.field('Right overhang', f"{vector['right']['overhang'] or '-'}")

    style.section('INSERT')
    style.field('Left enzyme', f"{insert['left']['enzyme']}")
    style.field('Left position', f"{insert['left']['position']}")
    style.field('Left recognition', f"{insert['left']['recognition_sequence']}")
    style.field('Left overhang', f"{insert['left']['overhang'] or '-'}")
    style.field('Right enzyme', f"{insert['right']['enzyme']}")
    style.field('Right position', f"{insert['right']['position']}")
    style.field('Right recognition', f"{insert['right']['recognition_sequence']}")
    style.field('Right overhang', f"{insert['right']['overhang'] or '-'}")

    style.section('CLONING JUNCTIONS')

    for side in (
        "left",
        "right",
    ):
        junction = analysis[
            "junctions"
        ][side]

        print(
            f"{side.upper():<7}"
            f"{junction['vector_enzyme']}"
            f" ({junction['vector_overhang']}) -> "
            f"{junction['insert_enzyme']}"
            f" ({junction['insert_overhang']}) "
            f"{junction['status']}"
        )
        print(
            f"        {junction['mechanism']}"
        )

    style.section('ORIENTATION')
    style.field('Intended orientation', f"{'COMPATIBLE' if validity['junctions_compatible'] else 'NOT VALID'}")
    style.field('Reverse orientation', f"{'COMPATIBLE' if validity['reverse_orientation_compatible'] else 'BLOCKED'}")
    style.field('Directional', f"{_yes_no(strategy['directional'])}")
    style.field('Self-ligation risk', f"{'POSSIBLE' if strategy['self_ligation_risk'] else 'REDUCED'}")

    style.section('DESIGN SCORE')
    style.field('Score', f"{score['score']}/{score['maximum']}")
    style.field('Category', f"{score['category']}")

    for component, value in score[
        "components"
    ].items():
        print(
            f"{component.replace('_', ' ').title():30}"
            f"{value:>3} pts"
        )

    style.section('VALIDATION')
    style.field('Required sites', f"{'PASS' if validity['required_sites'] else 'FAIL'}")
    style.field('Vector architecture', f"{'PASS' if validity['vector_unique'] else 'FAIL'}")
    style.field('Insert architecture', f"{'PASS' if validity['insert_unique'] else 'FAIL'}")
    style.field('Junctions', f"{'PASS' if validity['junctions_compatible'] else 'FAIL'}")
    style.field('No internal sites', f"{'PASS' if validity['no_internal_sites'] else 'FAIL'}")
    style.field('Overall valid', f"{'YES' if validity['valid'] else 'NO'}")

    print_design_reasons(
        analysis
    )

    ascii_orientation_visualization(
        analysis
    )

    style.section('INTERPRETATION')

    if validity[
        "valid"
    ]:
        print(
            "The selected design passes the implemented "
            "sequence-level cloning checks."
        )
    else:
        print(
            "The selected design does not pass all implemented "
            "sequence-level cloning checks."
        )

    print(
        "\nThis is a sequence-level design assessment. "
        "It does not guarantee experimental cleavage, ligation, "
        "enzyme activity, methylation behavior, buffer compatibility, "
        "or wet-lab success."
    )

    print(
        major_rule()
    )


def print_planner_report(
    planner,
):
    style.header('AUTOMATIC CLONING PLANNER')

    if not planner:
        print()
        style.note("No planner result available.")
        print(
            major_rule()
        )
        return

    print()
    style.field('Status', f"{planner.get('status', 'UNKNOWN')}")
    style.field('Strategies evaluated', f"{planner.get('evaluated_count', 0)}")
    style.field('Valid strategies', f"{planner.get('valid_count', 0)}")
    style.field('Rejected strategies', f"{planner.get('rejected_count', 0)}")

    best = planner.get(
        "best"
    )

    if best is None:
        diagnostics = planner.get(
            "diagnostics",
            {},
        )

        style.section('NO STRATEGY AVAILABLE')
        print(
            diagnostics.get(
                "primary_reason",
                "No cloning strategy could be evaluated.",
            )
        )

        recommendations = diagnostics.get(
            "recommendations",
            [],
        )

        if recommendations:
            print(
                "\nRecommendations:"
            )
            for recommendation in recommendations:
                print(
                    f"- {recommendation}"
                )

        for role in (
            "vector",
            "insert",
        ):
            section = diagnostics.get(
                role,
                {},
            )
            rows = section.get(
                "enzyme_diagnostics",
                [],
            )

            if not rows:
                continue

            print(
                f"\n{role.upper()} SITE DIAGNOSTICS"
            )
            print(
                minor_rule()
            )

            for row in rows:
                print(
                    f"{row['enzyme']:<12}"
                    f"{row['site_count']:>3} site(s)  "
                    f"{row['status']:<10}"
                    f"{row['sites']}"
                )

        print(
            major_rule()
        )
        return

    style.section('BEST STRATEGY')

    style.field('Vector', f"{best['vector']['left']['enzyme']} + "
        f"{best['vector']['right']['enzyme']}")
    style.field('Insert', f"{best['insert']['left']['enzyme']} + "
        f"{best['insert']['right']['enzyme']}")
    style.field('Score', f"{best['score']['score']}/"
        f"{best['score']['maximum']}")
    style.field('Category', f"{best['score']['category']}")
    style.field('Orientation', f"{best['strategy']['orientation']}")
    style.field('Directional', f"{_yes_no(best['strategy']['directional'])}")
    style.field('Valid', f"{_yes_no(best['validity']['valid'])}")

    style.section('JUNCTIONS')

    for side in (
        "left",
        "right",
    ):
        junction = best[
            "junctions"
        ][side]

        print(
            f"{side.upper():<7}"
            f"{junction['vector_enzyme']}"
            f" ({junction['vector_overhang']}) -> "
            f"{junction['insert_enzyme']}"
            f" ({junction['insert_overhang']}) "
            f"{junction['status']}"
        )

    style.section('ORIENTATION CHECK')

    style.field('Intended orientation', f"{'COMPATIBLE' if best['validity']['junctions_compatible'] else 'NOT VALID'}")
    style.field('Reverse orientation', f"{'COMPATIBLE' if best['validity']['reverse_orientation_compatible'] else 'BLOCKED'}")

    style.section('VALIDATION')

    validity = best[
        "validity"
    ]

    style.field('Required sites', f"{'PASS' if validity['required_sites'] else 'FAIL'}")
    style.field('Vector architecture', f"{'PASS' if validity['vector_unique'] else 'FAIL'}")
    style.field('Insert architecture', f"{'PASS' if validity['insert_unique'] else 'FAIL'}")
    style.field('Junctions', f"{'PASS' if validity['junctions_compatible'] else 'FAIL'}")
    style.field('No internal sites', f"{'PASS' if validity['no_internal_sites'] else 'FAIL'}")
    style.field('Overall valid', f"{'YES' if validity['valid'] else 'NO'}")

    print_design_reasons(
        best
    )

    style.section('RANKED STRATEGIES')

    for index, candidate in enumerate(
        planner.get(
            "strategies",
            [],
        ),
        start=1,
    ):
        vector = candidate[
            "vector"
        ]
        insert = candidate[
            "insert"
        ]
        strategy = candidate[
            "strategy"
        ]
        score = candidate[
            "score"
        ]

        vector_pair = (
            f"{vector['left']['enzyme']} + "
            f"{vector['right']['enzyme']}"
        )
        insert_pair = (
            f"{insert['left']['enzyme']} + "
            f"{insert['right']['enzyme']}"
        )

        print(
            f"{index:>2}. "
            f"V: {vector_pair:<22} "
            f"I: {insert_pair:<22} "
            f"{score['score']:>3}/100 "
            f"{'DIR' if strategy['directional'] else 'NON-DIR':<7} "
            f"{'VALID' if candidate['validity']['valid'] else 'REJECTED'}"
        )

    print(
        "\nImportant:"
    )
    print(
        "Ranking is based on sequence-level design properties only."
    )

    print(
        major_rule()
    )


print_advanced_cloning_report = (
    print_cloning_report
)


def _cloning_ui_select_insert(
    sequences,
    vector_index=None,
):
    if not sequences:
        print()
        style.note("No additional sequences are available.")
        return None, None, None

    records = []

    if isinstance(sequences, dict):
        records = list(sequences.items())
    else:
        for index, record in enumerate(sequences):
            if isinstance(record, dict):
                records.append((index, record))

    available = []

    style.header("SELECT INSERT")

    for key, record in records:
        sequence = record.get("sequence", "")
        header = record.get("header", str(key))

        if not sequence:
            continue

        if vector_index is not None and key == vector_index:
            continue

        available.append((key, header, sequence))

    if not available:
        style.note("No suitable insert sequence is available.")
        return None, None, None

    style.menu([("Insert", [(index, header) for index, (_, header, _) in
                            enumerate(available, start=1)])])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None, None, None

        for display_index, (key, header, sequence) in enumerate(
            available,
            start=1,
        ):
            if choice == str(display_index):
                return key, header, sequence

        style.invalid_selection()


def _cloning_ui_record_topology(
    sequences,
    key,
):
    """Topology of the record chosen as insert ("linear" if not recorded)."""
    try:
        record = sequences[key]
    except (KeyError, IndexError, TypeError):
        return "linear"

    if not isinstance(record, dict):
        return "linear"

    return record.get("topology", "linear")


def _cloning_ui_select_enzyme(
    prompt,
    allow_back=True,
):
    enzyme_names = list(RESTRICTION_ENZYMES.keys())

    style.header("SELECT RESTRICTION ENZYME")

    style.menu([("Enzymes", list(enumerate(enzyme_names, start=1)))],
               zero=style.BACK_LABEL if allow_back else None)

    while True:
        choice = style.ask_choice(prompt)

        if allow_back and choice == "0":
            return None

        try:
            index = int(choice) - 1
        except ValueError:
            style.invalid_selection()
            continue

        if 0 <= index < len(enzyme_names):
            return enzyme_names[index]

        style.invalid_selection()


def _cloning_ui_select_two_distinct_enzymes():
    first = _cloning_ui_select_enzyme(
        "Left enzyme"
    )

    if first is None:
        return None, None

    while True:
        second = _cloning_ui_select_enzyme(
            "Right enzyme"
        )

        if second is None:
            return None, None

        if second == first:
            print()
            style.error("The two enzymes must be different.")
            continue

        return first, second


def _cloning_ui_store_result(
    pipeline,
    result,
):
    if pipeline is None:
        return

    try:
        pipeline.context.store_cloning_analysis(
            result
        )
    except AttributeError:
        pass


def _cloning_ui_print_sequence_context(
    vector_header,
    vector_sequence,
    insert_header,
    insert_sequence,
    vector_topology="linear",
    insert_topology="linear",
):
    style.field('Vector', f"{vector_header}")
    style.field('Vector length', f"{len(normalize_sequence(vector_sequence))} bp ({vector_topology})")
    style.field('Insert', f"{insert_header}")
    style.field('Insert length', f"{len(normalize_sequence(insert_sequence))} bp ({insert_topology})")


def _cloning_ui_automatic_planner(
    vector_sequence,
    insert_sequence,
    pipeline=None,
    vector_topology="linear",
    insert_topology="linear",
):
    style.header(
        "AUTOMATIC CLONING PLANNER"
    )

    print(
        "\nSynapse evaluates the available enzyme architectures "
        "and ranks the resulting sequence-level strategies."
    )

    style.menu([
        (None, [('1', 'Find best strategy'), ('2', 'Show top strategies')]),
    ])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice in {"1", "2"}:
            break

        style.invalid_selection()

    if choice == "1":
        top_n = 1
    else:
        while True:
            raw = style.ask('Number of strategies to display [default 5]')

            if not raw:
                top_n = 5
                break

            try:
                top_n = int(raw)
                if top_n < 1:
                    raise ValueError
                break
            except ValueError:
                style.error("Please enter a positive integer.")

    result = automatic_cloning_planner(
        vector_sequence,
        insert_sequence,
        top_n=top_n,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )

    _cloning_ui_store_result(
        pipeline,
        result,
    )

    print_planner_report(result)
    style.pause()
    return result


def _cloning_ui_single_enzyme(
    vector_sequence,
    insert_sequence,
    pipeline=None,
    vector_topology="linear",
    insert_topology="linear",
):
    style.header(
        "SINGLE-ENZYME ANALYSIS"
    )

    print(
        "\nThe same restriction enzyme is evaluated "
        "for vector and insert."
    )

    style.menu([
        (None, [('1', 'Select enzyme')]),
    ])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            break

        style.invalid_selection()

    enzyme = _cloning_ui_select_enzyme(
        "Selection"
    )

    if enzyme is None:
        return None

    result = analyze_cloning(
        vector_sequence,
        insert_sequence,
        enzyme,
        enzyme,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )

    _cloning_ui_store_result(
        pipeline,
        result,
    )

    print_cloning_report(result)
    style.pause()
    return result


def _cloning_ui_select_two_enzyme_architecture(
    molecule_name,
):
    style.header(
        f"SELECT {molecule_name.upper()} ENZYMES"
    )

    print(
        f"\nSelect the two restriction enzymes for the "
        f"{molecule_name.lower()} ends."
    )
    print(
        "The two enzymes must be different."
    )
    print()

    left_enzyme = _cloning_ui_select_enzyme(
        f"{molecule_name} left enzyme"
    )

    if left_enzyme is None:
        return None, None

    right_enzyme = _cloning_ui_select_enzyme(
        f"{molecule_name} right enzyme"
    )

    if right_enzyme is None:
        return None, None

    if right_enzyme == left_enzyme:
        print()
        style.error("The two enzymes must be different.")
        style.pause()
        return _cloning_ui_select_two_enzyme_architecture(
            molecule_name
        )

    return left_enzyme, right_enzyme


def _cloning_ui_two_enzyme_strategy(
    vector_sequence,
    insert_sequence,
    pipeline=None,
    vector_topology="linear",
    insert_topology="linear",
):
    style.header(
        "TWO-ENZYME STRATEGY ANALYSIS"
    )

    print(
        "\nVector and insert are configured independently."
    )
    print(
        "This allows Synapse to evaluate four enzyme positions:"
    )
    print(
        "  Vector left / Vector right / Insert left / Insert right"
    )

    style.menu([
        (None, [('1', 'Configure enzyme architecture')]),
    ])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            break

        style.invalid_selection()

    vector_left, vector_right = (
        _cloning_ui_select_two_enzyme_architecture(
            "Vector"
        )
    )

    if vector_left is None:
        return None

    insert_left, insert_right = (
        _cloning_ui_select_two_enzyme_architecture(
            "Insert"
        )
    )

    if insert_left is None:
        return None

    style.header(
        "TWO-ENZYME STRATEGY"
    )

    print()
    style.field('Vector', f"{vector_left} + {vector_right}")
    style.field('Insert', f"{insert_left} + {insert_right}")

    style.menu([
        (None, [('1', 'Run analysis')]),
    ])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            break

        style.invalid_selection()

    result = evaluate_cloning_strategy(
        vector_sequence,
        insert_sequence,
        vector_left,
        vector_right,
        insert_left,
        insert_right,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )

    _cloning_ui_store_result(
        pipeline,
        result,
    )

    print_cloning_report(result)
    style.pause()
    return result


def run_cloning_analysis_ui(
    vector_sequence,
    insert_sequence=None,
    vector_header="Vector",
    insert_header="Insert",
    pipeline=None,
    sequences=None,
    vector_index=None,
    vector_topology="linear",
    insert_topology=None,
):
    """
    Interactive CLI for the canonical cloning engine.

    UI structure:
        CLONING ANALYSIS
            [1] Automatic Cloning Planner
            [2] Single-Enzyme Analysis
            [3] Two-Enzyme Strategy Analysis
            [0] Back

    Sequence-specific input and enzyme selection are requested only
    after the corresponding analysis is selected.

    The cloning engine itself is not modified. An insert chosen from
    ``sequences`` keeps the topology of its record; otherwise
    insert_topology (default linear) applies.
    """
    vector_validation = validate_dna_sequence(
        vector_sequence,
        "vector",
    )

    if not vector_validation["valid"]:
        style.header(
            "CLONING ANALYSIS"
        )
        print()
        style.error("Invalid vector sequence.")

        if vector_validation["invalid_positions"]:
            style.field('Invalid positions', f"{vector_validation['invalid_positions']}")

        print(minor_rule())
        return None

    if insert_sequence is None:
        insert_key, insert_header, insert_sequence = (
            _cloning_ui_select_insert(
                sequences,
                vector_index=vector_index,
            )
        )

        if insert_sequence is None:
            return None

        if insert_topology is None:
            insert_topology = _cloning_ui_record_topology(
                sequences,
                insert_key,
            )

    vector_topology = _normalize_topology(
        vector_topology
    )
    insert_topology = _normalize_topology(
        insert_topology or "linear"
    )

    insert_validation = validate_dna_sequence(
        insert_sequence,
        "insert",
    )

    if not insert_validation["valid"]:
        style.header(
            "CLONING ANALYSIS"
        )
        print()
        style.error("Invalid insert sequence.")

        if insert_validation["invalid_positions"]:
            style.field('Invalid positions', f"{insert_validation['invalid_positions']}")

        print(minor_rule())
        return None

    while True:
        style.header(
            "CLONING ANALYSIS"
        )

        _cloning_ui_print_sequence_context(
            vector_header,
            vector_sequence,
            insert_header,
            insert_sequence,
            vector_topology=vector_topology,
            insert_topology=insert_topology,
        )

        style.menu([
            ("Planning", [('1', 'Automatic Cloning Planner')]),
            ("Analysis", [('2', 'Single-Enzyme Analysis'), ('3', 'Two-Enzyme Strategy Analysis')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            result = _cloning_ui_automatic_planner(
                vector_sequence,
                insert_sequence,
                pipeline=pipeline,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )
            if result is not None:
                continue
            return None

        if choice == "2":
            result = _cloning_ui_single_enzyme(
                vector_sequence,
                insert_sequence,
                pipeline=pipeline,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )
            if result is not None:
                continue
            return None

        if choice == "3":
            result = _cloning_ui_two_enzyme_strategy(
                vector_sequence,
                insert_sequence,
                pipeline=pipeline,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )
            if result is not None:
                continue
            return None

        style.invalid_selection()


def run_cloning_ui_self_test():
    """The cloning dialog and its report printers are present (moved from analyses.cloning in P2-2)."""
    return all(
        callable(globals().get(name))
        for name in ("run_cloning_analysis_ui", "print_cloning_report", "print_planner_report")
    )
