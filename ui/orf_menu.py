"""Terminal dialogs of analyses/orf.py (moved here in P2-2).

The scientific functions stay in ``analyses.orf``; this module draws the
screens, reads the answers and calls them.
"""

from analyses.orf import (  # noqa: F401
    DEFAULT_MIN_ORF_LENGTH,
    DEFAULT_RANKING,
    DEFAULT_REPORT_NESTED,
    DEFAULT_START_CODON_POLICY,
    LEGACY_ORF_SETTINGS,
    _genetic_code_name,
    analyze_orf_relationships,
    analyze_orfs_enhanced,
    compare_orf_candidates,
    create_orf_preview,
    orf_end_on_sequence,
    summarize_orfs,
)
from ui import style
from ui.style import RULE_WIDTH, major_rule, minor_rule


def print_orf_report(
    orfs,
    min_length=DEFAULT_MIN_ORF_LENGTH
):

    style.header('ORF ANALYSIS')

    print()
    style.field('Minimum ORF Length', f"{min_length} bp")

    if not orfs:

        print("\nNo ORFs found.")

        style.section('ORF Definition')

        style.field('Search Space', "Both strands / six reading frames")

        style.field('Start Codons', "ATG, GTG, TTG")

        style.field('Stop Codons', "TAA, TAG, TGA")

        style.field('Minimum Length', f"{min_length} bp")

        print()
        print(major_rule())

        return

    # Infer sequence length from genomic coordinates.
    sequence_length = max(
        orf["stop"]
        for orf in orfs
    ) + 1

    summary = summarize_orfs(
        orfs,
        sequence_length
    )

    # ========================================================
    # SUMMARY
    # ========================================================

    style.section('ORF Summary')

    style.field('ORFs Found', f"{summary['total']}")

    style.field('Forward Strand', f"{summary['forward']}")

    style.field('Reverse Strand', f"{summary['reverse']}")

    style.field('Complete ORFs', f"{summary['complete']}")

    style.field('Partial ORFs', f"{summary['partial']}")

    style.field('Canonical Start (ATG)', f"{summary['canonical_start']}")

    style.field('Alternative Start', f"{summary['alternative_start']}")

    style.field('ORFs with Ambiguity', f"{summary['ambiguous']}")

    style.field('Translation-Safe', f"{summary['translation_safe']}")

    style.field('Kozak Context Available', f"{summary['kozak_available']}")

    style.field('Shine-Dalgarno Detected', f"{summary['shine_dalgarno_detected']}")

    style.field('Average ORF Length', f"{summary['average_length']:.1f} bp")

    style.field('Longest ORF', f"{summary['longest']['length']} bp")

    style.field('ORF Density', f"{summary['orf_density']:.2f}%")

    # ========================================================
    # FRAME STATISTICS
    # ========================================================

    style.section('Reading Frame Statistics')

    for frame in (
        "+1",
        "+2",
        "+3",
        "-1",
        "-2",
        "-3"
    ):

        statistics = (
            summary["frame_statistics"]
            [frame]
        )

        print(
            f"Frame {frame}: "
            f"{statistics['count']} ORFs | "
            f"{statistics['complete']} complete | "
            f"{statistics['partial']} partial"
        )

        if statistics["longest"]:

            print(
                f"           Longest: "
                f"{statistics['longest']['length']} bp"
            )

    # ========================================================
    # SEARCH PARAMETERS
    # ========================================================

    style.section('Search Parameters')

    style.field('Search Space', "Both strands / six reading frames")

    style.field('Start Codons', "ATG, GTG, TTG")

    style.field('Canonical Start', "ATG")

    style.field('Alternative Starts', "GTG, TTG")

    style.field('Stop Codons', "TAA, TAG, TGA")

    style.field('Minimum Length', f"{min_length} bp")

    style.field('Stop Rule', "First definite downstream "
        "in-frame stop")

    style.field('Ambiguous Starts', "Only definite starts initiate ORFs")

    style.field('Ambiguous Stops', "Possible-only stops do not terminate ORFs")

    # ========================================================
    # RANKED ORFS
    # ========================================================

    style.section('Ranked ORF Candidates')

    for orf in orfs:

        print()

        print(
            f"ORF #{orf['rank']}"
        )

        print(minor_rule())

        style.field('Heuristic score', f"{orf['heuristic_score']['score']}/100 "
            f"({orf['heuristic_score']['classification']})")

        style.field('CDS Quality', f"{orf['cds_quality']['score']}/100 "
            f"({orf['cds_quality']['classification']})")

        style.field('Strand', f"{orf['strand']}")

        style.field('Reading Frame', f"{orf['frame']:+d}")

        style.field('Start Position', f"{orf['start'] + 1}")

        style.field('Stop Position', f"{orf_end_on_sequence(orf)}")

        if orf.get("wraps_origin"):

            style.field('Across origin', 'yes (circular sequence)')

        style.field('Length', f"{orf['length']} bp")

        style.field('Codons', f"{orf['codon_count']}")

        style.field('Potential Protein', f"{orf['protein_length']} aa")

        style.field('Start Codon', f"{orf['start_codon']}")

        style.field('Start Type', f"{orf['start_type']}")

        style.field('Status', f"{orf['status']}")

        if orf["stop_codon"] is None:

            style.field('Stop Codon', "Not detected")

        else:

            style.field('Stop Codon', f"{orf['stop_codon']}")

        # ====================================================
        # START CONTEXT
        # ====================================================

        style.section('Start Context')

        kozak = (
            orf["start_context"]
            ["kozak"]
        )

        if kozak["available"]:

            style.field('Kozak Context', f"{kozak['context']}")

            style.field('Kozak Score', f"{kozak['score']}/3")

            print(
                f"Kozak Classification:"
                f" {kozak['classification']}"
            )

        else:

            style.field('Kozak Context', "Unavailable")

        shine = (
            orf["start_context"]
            ["shine_dalgarno"]
        )

        if shine["detected"]:

            style.field('Shine-Dalgarno', f"{shine['motif']}")

            style.field('SD Distance', f"{shine['distance_to_start']} nt")

        else:

            style.field('Shine-Dalgarno', "Not detected")

        # ====================================================
        # CODON STATISTICS
        # ====================================================

        codon_usage = (
            orf["codon_usage"]
        )

        style.section('Codon Statistics')

        style.field('Total Codons', f"{codon_usage['total_codons']}")

        style.field('Unique Codons', f"{codon_usage['unique_codons']}")

        style.field('Codon Diversity', f"{orf['codon_diversity']:.3f}")

        if (
            codon_usage["gc3_content"]
            is not None
        ):

            style.field('GC3 Content', f"{codon_usage['gc3_content']:.2f}%")

        else:

            style.field('GC3 Content', "Not available")

        # ====================================================
        # AMBIGUITY
        # ====================================================

        style.section('Sequence Ambiguity')

        style.field('Ambiguous Codons', f"{orf['ambiguous_codon_count']}")

        style.field('Deterministic', f"{orf['deterministic_ambiguous_count']}")

        print(
            f"Translation-Uncertain:"
            f" {orf['uncertain_codon_count']}"
        )

        print(
            f"Possible Stop Codons:"
            f" {orf['possible_stop_codon_count']}"
        )

        if orf["translation_safe"]:

            if (
                orf["ambiguous_codon_count"]
                == 0
            ):

                style.field('Translation Status', " Unambiguous")

            else:

                style.field('Translation Status', " Deterministically resolvable")

        else:

            style.field('Translation Status', " Contains uncertainty")

        # ====================================================
        # AMBIGUOUS CODON DETAILS
        # ====================================================

        if orf["ambiguous_codons"]:

            style.section('Ambiguous Codon Details')

            for item in (
                orf["ambiguous_codons"]
            ):

                codon_number = (
                    (
                        item["position"] - 1
                    )
                    // 3
                ) + 1

                if (
                    item[
                        "translation_status"
                    ]
                    == "Deterministic"
                ):

                    print(
                        f"Codon {codon_number}: "
                        f"{item['codon']} -> "
                        f"{item['translation']}"
                    )

                else:

                    possibilities = (
                        "/".join(
                            item[
                                "possibilities"
                            ]
                        )
                    )

                    if item[
                        "contains_stop"
                    ]:

                        possibilities += (
                            "/STOP"
                        )

                    print(
                        f"Codon {codon_number}: "
                        f"{item['codon']} -> "
                        f"{possibilities}"
                    )

        # ====================================================
        # SEQUENCE
        # ====================================================

        style.section('Sequence')

        print(
            create_orf_preview(
                orf["sequence"]
            )
        )

    # ========================================================
    # INTERPRETATION
    # ========================================================

    style.section('Interpretation')

    print(
        "ORFs are sequence-defined candidate "
        "coding regions."
    )

    print(
        "ORF detection alone does not prove "
        "biological expression or protein function."
    )

    print(
        "ATG is treated as the canonical start codon, "
        "while GTG and TTG are treated as alternative "
        "start codons."
    )

    print(
        "The analysis evaluates all six reading frames "
        "on both DNA strands."
    )

    print(
        "IUPAC ambiguity is explicitly evaluated rather "
        "than silently ignored."
    )

    print(
        "Kozak and Shine-Dalgarno results are sequence-based "
        "context indicators and do not prove translation."
    )

    print(
        "Heuristic and CDS quality scores are "
        "ranking systems, not experimental evidence."
    )

    print()

    print(major_rule())


def print_orf_candidate_comparison(orfs, top_n=5):
    """Print a compact comparison of the highest-ranked candidates."""
    candidates = compare_orf_candidates(orfs, top_n=top_n)

    style.section('Top ORF Candidate Comparison')

    if not candidates:
        style.note("No ORF candidates available.")
        return

    for candidate in candidates:
        print(
            f"#{candidate['rank']} | "
            f"Frame {candidate['frame']:+d} | "
            f"{candidate['length']} bp | "
            f"{candidate['protein_length']} aa | "
            f"{candidate['start_codon']} | "
            f"{candidate['status']} | "
            f"Heuristic {candidate['heuristic_score']}/100 | "
            f"CDS {candidate['cds_quality_score']}/100"
        )


def print_orf_relationship_report(orfs):
    """Print detected overlap/nesting relationships."""
    relationships = analyze_orf_relationships(orfs)

    style.section('ORF Relationships')

    style.field('Overlapping pairs', f"{relationships['overlap_count']}")
    style.field('Nested pairs', f"{relationships['nested_count']}")

    if not relationships["relationships"]:
        print("No overlapping or nested ORF relationships detected.")
        return

    for item in relationships["relationships"]:
        print(
            f"ORF #{item['orf_a_rank']} ↔ ORF #{item['orf_b_rank']} | "
            f"{item['relation']} | "
            f"overlap {item['overlap_length']} bp | "
            f"same strand: {'YES' if item['same_strand'] else 'NO'} | "
            f"same frame: {'YES' if item['same_frame'] else 'NO'}"
        )


def _orf_status_label(orf):
    return "COMPLETE" if orf["status"] == "Complete" else "PARTIAL"


def print_orf_summary_card(analysis):
    """Compact executive summary for the ORF UI."""
    style.header("ORF ANALYSIS SUMMARY")

    summary = analysis.get("summary", {})
    top = analysis.get("top_candidate")

    style.field('Sequence length', f"{analysis.get('sequence_length', 0)} bp")
    style.field('ORF candidates', f"{summary.get('total', 0)}")

    if top:
        style.field('Best candidate', f"ORF #{top['rank']}")
        style.field('Heuristic score', f"{top['heuristic_score']['score']}/100 "
            f"({top['heuristic_score']['classification']})")
        style.field('CDS quality', f"{top['cds_quality']['score']}/100 "
            f"({top['cds_quality']['classification']})")
        style.field('Status', f"{top['status']}")
    else:
        style.field('Best candidate', "None")

    print()
    style.field('Candidates', f"{summary.get('complete', 0)} complete, "
        f"{summary.get('partial', 0)} partial; "
        f"{summary.get('forward', 0)} forward, "
        f"{summary.get('reverse', 0)} reverse")


def print_orf_candidates(orfs, top_n=None):
    """Display a clean candidate list without asking for input."""
    style.header("ORF CANDIDATES")

    if not orfs:
        style.note("No ORF candidates available.")
        return

    candidates = orfs if top_n is None else orfs[:max(0, top_n)]

    widths = (8, 6, 9, 8, 6, 10)
    print(style.columns(
        ("ID", "Frame", "Length", "Protein", "Start", "Status", "Heuristic"), widths))
    print(minor_rule())

    for orf in candidates:
        print(style.columns((
            f"ORF-{orf['rank']:03d}",
            f"{orf['frame']:+d}",
            f"{orf['length']} bp",
            f"{orf['protein_length']} aa",
            orf["start_codon"],
            _orf_status_label(orf),
            f"{orf['heuristic_score']['score']}/100",
        ), widths))


def print_orf_candidate_details(orf):
    """Detailed view for one candidate."""
    style.header(f"ORF-{orf['rank']:03d} — CANDIDATE DETAILS")

    style.field('Identifier', f"ORF-{orf['rank']:03d}")
    style.field('Rank', f"#{orf['rank']}")
    style.field('Strand', f"{orf['strand']}")
    style.field('Reading frame', f"{orf['frame']:+d}")
    style.field('Genomic start', f"{orf['start'] + 1}")
    style.field('Genomic stop', f"{orf_end_on_sequence(orf)}")
    if orf.get("wraps_origin"):
        style.field('Across origin', 'yes (circular sequence)')
    style.field('Length', f"{orf['length']} bp")
    style.field('Potential protein', f"{orf['protein_length']} aa")
    style.field('Start codon', f"{orf['start_codon']}")
    style.field('Start type', f"{orf['start_type']}")
    style.field('Status', f"{orf['status']}")

    style.section('Quality')
    style.field('Heuristic score', f"{orf['heuristic_score']['score']}/100 "
        f"({orf['heuristic_score']['classification']})")
    style.field('CDS quality', f"{orf['cds_quality']['score']}/100 "
        f"({orf['cds_quality']['classification']})")
    style.field('Translation safe', f"{'YES' if orf['translation_safe'] else 'NO'}")

    style.section('Start Context')

    kozak = orf["start_context"]["kozak"]
    if kozak["available"]:
        style.field('Kozak context', f"{kozak['context']}")
        style.field('Kozak score', f"{kozak['score']}/3")
        style.field('Kozak classification', f"{kozak['classification']}")
    else:
        style.field('Kozak context', "Unavailable")

    shine = orf["start_context"]["shine_dalgarno"]
    if shine["detected"]:
        style.field('Shine-Dalgarno', f"{shine['motif']}")
        style.field('SD distance', f"{shine['distance_to_start']} nt")
    else:
        style.field('Shine-Dalgarno', "Not detected")

    style.section('Sequence')
    print(create_orf_preview(orf["sequence"]))


def print_orf_comparison_table(orfs, top_n=5):
    """Display the highest-ranked candidates side by side."""
    candidates = compare_orf_candidates(orfs, top_n=top_n)

    style.header("COMPARE ORFs")

    if not candidates:
        style.note("No ORF candidates available.")
        return

    print(
        f"{'Feature':<24}"
        + "".join(f"ORF #{c['rank']:<14}" for c in candidates)
    )
    # A table rule spans its own table; the standard width is the floor.
    print(minor_rule() if 24 + 14 * len(candidates) <= RULE_WIDTH
          else "-" * (24 + 14 * len(candidates)))

    rows = [
        ("Strand", lambda c: c["strand"]),
        ("Reading frame", lambda c: f"{c['frame']:+d}"),
        ("Length", lambda c: f"{c['length']} bp"),
        ("Protein", lambda c: f"{c['protein_length']} aa"),
        ("Start codon", lambda c: c["start_codon"]),
        ("Status", lambda c: c["status"]),
        ("Heuristic score", lambda c: f"{c['heuristic_score']}/100"),
        ("CDS quality", lambda c: f"{c['cds_quality_score']}/100"),
        (
            "Translation safe",
            lambda c: "YES" if c["translation_safe"] else "NO",
        ),
        ("Ambiguous codons", lambda c: str(c["ambiguous_codons"])),
    ]

    for label, getter in rows:
        values = "".join(f"{getter(c):<14}" for c in candidates)
        print(f"{label:<24}{values}")


def print_orf_relationships(orfs):
    """Display sequence-coordinate relationships between ORFs."""
    data = analyze_orf_relationships(orfs)

    style.header("ORF RELATIONSHIPS")

    style.field('Relationships found', f"{data['relationship_count']}")
    style.field('Overlapping pairs', f"{data['overlap_count']}")
    style.field('Nested pairs', f"{data['nested_count']}")

    if not data["relationships"]:
        print("\nNo overlapping or nested ORF relationships detected.")
        return

    print()

    for item in data["relationships"]:
        print(
            f"ORF #{item['orf_a_rank']} <-> ORF #{item['orf_b_rank']}"
        )
        style.field('Relationship', f"{item['relation']}")
        style.field('Overlap', f"{item['overlap_length']} bp")
        style.field('ORF A coverage', f"{item['overlap_fraction_a'] * 100:.1f}%")
        style.field('ORF B coverage', f"{item['overlap_fraction_b'] * 100:.1f}%")
        style.field('Same strand', f"{'YES' if item['same_strand'] else 'NO'}")
        style.field('Same frame', f"{'YES' if item['same_frame'] else 'NO'}")
        print(minor_rule())


def select_orf(orfs):
    """
    Single source of truth for ORF selection.

    Important:
    The caller must not call another selection prompt after this function
    returns. The selected ORF is already the object used by downstream RNA
    and protein analysis.
    """
    if not orfs:
        print()
        style.note("No ORF candidates available for selection.")
        return None

    style.header("SELECT ORF")

    style.menu([("Available candidates", [
        (index,
         f"ORF-{orf['rank']:03d} | "
         f"Frame {orf['frame']:+d} | "
         f"{orf['length']} bp | "
         f"{orf['protein_length']} aa | "
         f"{orf['start_codon']} | "
         f"{orf['status']} | "
         f"Heuristic {orf['heuristic_score']['score']}/100")
        for index, orf in enumerate(orfs, start=1)
    ])])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice.isdigit():
            index = int(choice)

            if 1 <= index <= len(orfs):
                selected = orfs[index - 1]

                style.header("ORF SELECTED")

                style.field('Identifier', f"ORF-{selected['rank']:03d}")
                style.field('Rank', f"#{selected['rank']}")
                style.field('Strand', f"{selected['strand']}")
                style.field('Reading frame', f"{selected['frame']:+d}")
                style.field('Length', f"{selected['length']} bp")
                style.field('Protein length', f"{selected['protein_length']} aa")
                style.field('Start codon', f"{selected['start_codon']}")
                style.field('Status', f"{selected['status']}")
                style.field('Heuristic score', f"{selected['heuristic_score']['score']}/100")
                style.field('CDS quality', f"{selected['cds_quality']['score']}/100")

                style.section('Pipeline Status')
                style.success("ORF Selected")
                style.success("Ready for RNA Analysis")
                style.success("Ready for Protein Analysis")
                style.success("Ready for CRISPR guide design.")

                return selected

        style.invalid_selection()
        style.note("Enter a candidate number or 0 to go back.")


def _partial_label(settings):
    if not settings["include_partial"]:
        return "excluded"
    if settings.get("include_5prime_partial"):
        return "5' and 3' ends"
    return "3' end (missing stop)"


def _orf_settings_menu(settings):
    """Let the user change the ORF search policy. Returns updated settings."""
    settings = dict(settings)
    while True:
        style.header("ORF SEARCH SETTINGS")
        style.menu([
            ("Search policy", [
                (1, f"Minimum ORF length:  {settings['min_length']} nt (stop codon included)"),
                (2, f"Start codons:        {settings['start_codon_policy']}"),
                (3, f"Nested ORFs:         {'reported' if settings['nested'] else 'longest per stop only'}"),
                (4, f"Ranking:             {settings['ranking']}"),
                (5, f"Partial ORFs:        {_partial_label(settings)}"),
                (6, f"Genetic code:        {_genetic_code_name(settings.get('genetic_code', 1))}"),
            ]),
            ("Defaults", [
                (7, "Restore Synapse defaults (ORFfinder-style)"),
                (8, "Restore legacy settings (before 2026-10-29)"),
            ]),
        ])
        style.note("Defaults follow NCBI ORFfinder conventions (75 nt, ATG only, one ORF per "
                   "stop). They are screening conventions, not biological thresholds; adjust "
                   "them to your organism and question.")
        style.note("Synapse counts the stop codon in the minimum length, ORFfinder does not: "
                   "78 nt here gives ORFfinder's 75 nt.")

        choice = style.ask_choice()

        if choice == "0":
            return settings
        if choice == "1":
            raw = style.ask('Minimum length in nt [75]')
            try:
                value = int(raw) if raw else 75
                if value < 3:
                    raise ValueError
                settings["min_length"] = value
            except ValueError:
                style.error("Please enter an integer >= 3.")
        elif choice == "2":
            settings["start_codon_policy"] = (
                "atg_and_alternative"
                if settings["start_codon_policy"] == "atg_only"
                else "atg_only"
            )
        elif choice == "3":
            settings["nested"] = not settings["nested"]
        elif choice == "4":
            settings["ranking"] = (
                "heuristic" if settings["ranking"] == "length" else "length"
            )
        elif choice == "5":
            # 3' end only -> both ends -> excluded -> 3' end only
            if settings["include_partial"] and not settings.get("include_5prime_partial"):
                settings["include_5prime_partial"] = True
            elif settings["include_partial"]:
                settings["include_partial"] = False
                settings["include_5prime_partial"] = False
            else:
                settings["include_partial"] = True
        elif choice == "7":
            settings.update(
                min_length=DEFAULT_MIN_ORF_LENGTH,
                start_codon_policy=DEFAULT_START_CODON_POLICY,
                nested=DEFAULT_REPORT_NESTED,
                ranking=DEFAULT_RANKING,
                genetic_code=1,
                include_5prime_partial=False,
            )
        elif choice == "8":
            settings.update(LEGACY_ORF_SETTINGS)
            settings["genetic_code"] = 1
            settings["include_5prime_partial"] = False
        elif choice == "6":
            from core.seqtools import available_genetic_codes
            for table in available_genetic_codes():
                print(f"  {table.table_id:>2}  {table.name}")
            raw = style.ask('NCBI table number [1]')
            try:
                value = int(raw) if raw else 1
                _genetic_code_name(value)
                settings["genetic_code"] = value
            except ValueError as exc:
                style.error(f"Genetic code not available: {exc}")
        else:
            style.invalid_selection()


def run_orf_analysis_ui(
    dna,
    include_partial=True,
    min_length=DEFAULT_MIN_ORF_LENGTH,
    start_codon_policy=None,
    nested=DEFAULT_REPORT_NESTED,
    ranking=DEFAULT_RANKING,
    genetic_code=1,
    topology="linear",
):
    """
    Professional ORF module.

    Workflow:
        Settings -> Discovery -> Candidates -> Compare -> Relationships -> Select

    Returns the selected ORF or None. The selected ORF carries the search
    parameters that produced it under "search_parameters".
    """
    settings = {
        "include_partial": include_partial,
        "include_5prime_partial": False,
        "min_length": min_length,
        "start_codon_policy": start_codon_policy or DEFAULT_START_CODON_POLICY,
        "nested": nested,
        "ranking": ranking,
        "genetic_code": genetic_code,
    }

    def run_search():
        return analyze_orfs_enhanced(
            dna,
            include_partial=settings["include_partial"],
            min_length=settings["min_length"],
            start_codon_policy=settings["start_codon_policy"],
            nested=settings["nested"],
            ranking=settings["ranking"],
            genetic_code=settings["genetic_code"],
            topology=topology,
            include_5prime_partial=settings["include_5prime_partial"],
        )

    analysis = run_search()

    while True:
        style.header("ORF ANALYSIS")

        style.field('Sequence length', f"{analysis['sequence_length']} bp")
        if str(topology).lower() == "circular":
            style.field('Topology', 'circular (ORFs across the origin included)')
        style.field('ORF candidates', f"{analysis['candidate_count']}")
        style.field('Search', f">= {settings['min_length']} nt, "
            f"{settings['start_codon_policy']}, "
            f"{'nested' if settings['nested'] else 'longest per stop'}, "
            f"ranked by {settings['ranking']}")

        if settings["genetic_code"] != 1:
            style.field('Genetic code', f"{_genetic_code_name(settings['genetic_code'])}")

        top = analysis.get("top_candidate")

        if top:
            style.field('Best candidate', f"ORF-{top['rank']:03d} "
                f"({top['length']} bp, {top['start_codon']}, {top['status']})")
        else:
            style.field('Best candidate', "None")

        style.menu([
            ("Discovery", [(1, "ORF Discovery"), (2, "ORF Candidates")]),
            ("Evaluation", [(3, "Compare ORFs"), (4, "ORF Relationships")]),
            ("Selection", [(5, "Select ORF")]),
            ("Settings", [(6, "Search Settings")]),
        ])

        choice = style.ask_choice()

        if choice == "1":
            print_orf_report(
                analysis["orfs"],
                min_length=settings["min_length"],
            )
            style.pause()

        elif choice == "2":
            print_orf_summary_card(analysis)
            print_orf_candidates(analysis["orfs"])
            style.pause()

        elif choice == "3":
            print_orf_comparison_table(analysis["orfs"])
            style.pause()

        elif choice == "4":
            print_orf_relationships(analysis["orfs"])
            style.pause()

        elif choice == "5":
            selected = select_orf(analysis["orfs"])
            if selected is not None:
                selected["search_parameters"] = dict(analysis.get("parameters", {}))
                return selected

        elif choice == "6":
            settings = _orf_settings_menu(settings)
            analysis = run_search()

        elif choice == "0":
            return None

        else:
            style.invalid_selection()


show_orf_candidates = print_orf_candidates


show_orf_comparison = print_orf_comparison_table


show_orf_relationships = print_orf_relationships
