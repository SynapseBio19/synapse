"""Terminal dialogs of analyses/mutation.py (moved here in P2-2).

The scientific functions stay in ``analyses.mutation``; this module draws the
screens, reads the answers and calls them.
"""

from __future__ import annotations

from typing import Any
from typing import List
from typing import Tuple
import math
from analyses.mutation import (  # noqa: F401
    CODING_FRAME_ASSUMPTION,
    MutationAnalysisError,
    MutationAnalysisResult,
    SequencePromptCancelled,
    annotate_mutations_with_coding_consequences,
    assess_coding_frame,
    build_mutation_report,
    calculate_mutation_statistics,
    compare_sequences,
    normalize_dna,
    normalize_iupac_sequence,
    orf_reference_cds,
    run_mutation_self_test,
)
from ui import style
from ui.style import minor_rule


def print_mutation_report(
    result: MutationAnalysisResult,
    *,
    include_alignment: bool = False,
) -> None:

    print(
        build_mutation_report(
            result,
            include_alignment=include_alignment,
        )
    )


def _extract_sequence_record(
    record_key: Any,
    record: Any,
) -> Tuple[str, Any]:

    if isinstance(record, dict):

        header = str(
            record.get("header")
            or record.get("name")
            or record.get("id")
            or record_key
        )

        sequence = (
            record.get("sequence")
            or record.get("dna")
        )

        return header, sequence

    if (
        isinstance(record, (tuple, list))
        and len(record) >= 2
    ):

        return (
            str(record[0]),
            record[1],
        )

    return (
        str(record_key),
        record,
    )


def _iter_sequence_records(
    sequences: Any,
) -> List[Tuple[str, str]]:

    available = []

    if isinstance(
        sequences,
        dict,
    ):

        iterator = sequences.items()

    elif isinstance(
        sequences,
        (list, tuple),
    ):

        iterator = enumerate(
            sequences,
            1,
        )

    else:

        return available

    for key, record in iterator:

        try:

            header, sequence = (
                _extract_sequence_record(
                    key,
                    record,
                )
            )

            if sequence is None:
                continue

            sequence = str(
                sequence
            ).strip()

            if sequence:

                available.append(
                    (
                        str(header),
                        sequence,
                    )
                )

        except Exception:
            continue

    return available


def _prompt_sequence(
    prompt: str,
    sequences=None,
) -> Tuple[str, str]:

    available = _iter_sequence_records(
        sequences
    )

    if available:

        # [0] raises SequencePromptCancelled; both menu branches that ask
        # for sequences catch it and return to the mutation menu.
        style.menu([
            (f"Available sequences for {prompt}",
             [(i, header) for i, (header, _) in enumerate(available, 1)]),
            (None, [("M", "Enter sequence manually")]),
        ], zero=style.CANCEL_LABEL)

        while True:

            choice = style.ask_choice(f"{prompt} selection")

            if choice == "0":
                raise SequencePromptCancelled(prompt)

            if choice.upper() == "M":
                break

            if (
                choice.isdigit()
                and 1 <= int(choice)
                <= len(available)
            ):

                header, sequence = (
                    available[
                        int(choice) - 1
                    ]
                )

                return (
                    header,
                    normalize_iupac_sequence(
                        sequence
                    ),
                )

            style.invalid_selection()

    header = style.ask(f"{prompt} name")

    if not header:
        header = prompt

    sequence = style.ask(f"{prompt} DNA sequence")

    return (
        header,
        normalize_iupac_sequence(
            sequence
        ),
    )


def _store_pipeline_result(
    pipeline,
    key: str,
    value: Any,
) -> None:

    if pipeline is None:
        return

    try:

        if (
            hasattr(
                pipeline,
                "store",
            )
            and callable(
                pipeline.store
            )
        ):

            pipeline.store(
                key,
                value,
            )

        elif (
            hasattr(
                pipeline,
                "results",
            )
            and isinstance(
                pipeline.results,
                dict,
            )
        ):

            pipeline.results[
                key
            ] = value

        elif isinstance(
            pipeline,
            dict,
        ):

            pipeline[key] = value

        else:
            import logging
            logging.getLogger("synapse.mutation").warning(
                "Pipeline object has no result store; %s result not stored.", key
            )

    except Exception:
        # Storing must never destroy a finished analysis, but the failure
        # is logged instead of being hidden (audit Phase V).
        import logging
        logging.getLogger("synapse.mutation").exception(
            "Could not store %s result in the pipeline.", key
        )


def run_mutation_analysis_ui(
    dna=None,
    pipeline=None,
    header="Sequence",
    sequences=None,
):

    last_result = None

    while True:

        style.header('MUTATION ANALYSIS')
        style.field("Reference", header if dna else "none loaded")
        if dna:
            style.field("Length", f"{len(dna)} bp")

        style.menu([
            ("Comparison", [('1', 'Compare Sequences')]),
            ("Analysis", [('2', 'Mutation Statistics'), ('3', 'Coding Consequence Analysis')]),
            ("Output", [('4', 'Mutation Report')]),
            ("System", [('5', 'Run Mutation Self-Test')]),
        ])

        choice = style.ask_choice()

        # ------------------------------------------------------------------
        # BACK
        # ------------------------------------------------------------------

        if choice == "0":
            return last_result

        # ------------------------------------------------------------------
        # COMPARE
        # ------------------------------------------------------------------

        if choice == "1":

            try:

                if dna:

                    reference_name = (
                        header
                    )

                    reference_sequence = (
                        normalize_iupac_sequence(
                            dna
                        )
                    )

                else:

                    (
                        reference_name,
                        reference_sequence,
                    ) = _prompt_sequence(
                        "Reference",
                        sequences,
                    )

                (
                    variant_name,
                    variant_sequence,
                ) = _prompt_sequence(
                    "Variant",
                    sequences,
                )

                last_result = compare_sequences(
                    reference_sequence,
                    variant_sequence,
                    reference_name=reference_name,
                    variant_name=variant_name,
                )

                # Coding annotations are only meaningful when the user
                # supplies/uses CDS sequences. Do not silently claim that
                # an arbitrary DNA sequence is coding.
                last_result.statistics = (
                    calculate_mutation_statistics(
                        last_result.mutations,
                        reference_length=(
                            last_result.reference_length
                        ),
                    )
                )

                print()

                print_mutation_report(
                    last_result
                )

                _store_pipeline_result(
                    pipeline,
                    "mutation",
                    last_result.to_dict(),
                )

            except SequencePromptCancelled:

                continue

            except MutationAnalysisError as error:

                style.error(f"{error}")

            except Exception as error:

                style.error(f"Mutation analysis failed: {error}")

        # ------------------------------------------------------------------
        # STATISTICS
        # ------------------------------------------------------------------

        elif choice == "2":

            if last_result is None:

                style.note("Run a sequence comparison first.")

                continue

            s = last_result.statistics

            if (
                s.transition_transversion_ratio
                == math.inf
            ):

                ratio = "inf"

            elif (
                s.transition_transversion_ratio
                is None
            ):

                ratio = "N/A"

            else:

                ratio = (
                    f"{s.transition_transversion_ratio:.3f}"
                )

            style.header('MUTATION STATISTICS')

            style.field('Total mutations', f"{s.total_mutations}")

            style.field('Alignment identity', f"{last_result.identity_percent:.2f}%")

            style.field('Sequence similarity', f"{last_result.sequence_similarity_percent:.2f}%")

            style.field('SNPs', f"{s.snps}")

            style.field('Insertions', f"{s.insertions}")

            style.field('Deletions', f"{s.deletions}")

            style.field('Transitions', f"{s.transitions}")

            style.field('Transversions', f"{s.transversions}")

            style.field('Ti/Tv ratio', f"{ratio}")

            style.field('Net length change', f"{s.net_length_change:+d} bp")

            style.field('Mutation density', f"{s.mutation_density_per_kb:.2f}/kb")

            style.field('Synonymous', f"{s.synonymous}")

            style.field('Missense', f"{s.missense}")

            style.field('Nonsense', f"{s.nonsense}")

            style.field('Frameshifts', f"{s.frameshifts}")

            style.field('In-frame indels', f"{s.in_frame_indels}")

        # ------------------------------------------------------------------
        # CODING CONSEQUENCE
        # ------------------------------------------------------------------

        elif choice == "3":

            try:

                if last_result is None:

                    if dna:

                        ref_name = header

                        ref = normalize_dna(
                            dna,
                            allow_n=False,
                        )

                    else:

                        (
                            ref_name,
                            ref,
                        ) = _prompt_sequence(
                            "Reference",
                            sequences,
                        )

                    (
                        var_name,
                        var,
                    ) = _prompt_sequence(
                        "Variant",
                        sequences,
                    )

                    last_result = compare_sequences(
                        ref,
                        var,
                        reference_name=ref_name,
                        variant_name=var_name,
                    )

                orf_ref_cds = orf_var_cds = None
                selected_orf = getattr(
                    getattr(pipeline, "context", None), "selected_orf", None
                )
                if selected_orf and dna:
                    if style.confirm("Use the selected pipeline ORF as CDS?", default_yes=True):
                        reference_plain = last_result.aligned_reference.replace("-", "")
                        variant_plain = last_result.aligned_variant.replace("-", "")
                        orf_ref_cds, orf_var_cds, orf_text = orf_reference_cds(
                            selected_orf, reference_plain, variant_plain
                        )
                        print(f"CDS from {orf_text}")
                        if orf_var_cds is None:
                            print(
                                "Variant length differs from the reference: enter the "
                                "variant CDS explicitly (ORF coordinates do not carry over)."
                            )
                ref_cds = orf_ref_cds or style.ask('Reference CDS (press Enter to use reference sequence)')

                if not ref_cds:

                    ref_cds = (
                        last_result
                        .aligned_reference
                        .replace(
                            "-",
                            "",
                        )
                    )

                var_cds = orf_var_cds or style.ask('Variant CDS (press Enter to use variant sequence)')

                if not var_cds:

                    var_cds = (
                        last_result
                        .aligned_variant
                        .replace(
                            "-",
                            "",
                        )
                    )

                coding = (
                    annotate_mutations_with_coding_consequences(
                        last_result.mutations,
                        ref_cds,
                        var_cds,
                    )
                )

                frame_warnings = assess_coding_frame(ref_cds)
                coding.notes.append(CODING_FRAME_ASSUMPTION)
                coding.notes.extend(
                    f"Reference CDS: {warning}" for warning in frame_warnings
                )
                if frame_warnings:
                    print()
                    for warning in frame_warnings:
                        style.warning(f"Coding frame assumption: {warning}")
                last_result.coding_analysis = (
                    coding
                )

                last_result.statistics = (
                    calculate_mutation_statistics(
                        last_result.mutations,
                        reference_length=(
                            last_result.reference_length
                        ),
                    )
                )

                style.header('CODING CONSEQUENCE ANALYSIS')

                style.field('Reference protein', f"{coding.reference_protein}")

                style.field('Variant protein', f"{coding.alternative_protein}")

                style.field('Synonymous', f"{coding.synonymous_count}")

                style.field('Missense', f"{coding.missense_count}")

                style.field('Nonsense', f"{coding.nonsense_count}")

                style.field('Start-loss', f"{coding.start_loss_count}")

                style.field('Stop-loss', f"{coding.stop_loss_count}")

                style.field('Frameshift', f"{'YES' if coding.frameshift else 'NO'}")

                if (
                    coding.frameshift_position
                    is not None
                ):

                    style.field('Frameshift position', f"{coding.frameshift_position}")

                if (
                    coding.premature_stop_position
                    is not None
                ):

                    style.field('Premature stop', f"codon "
                        f"{coding.premature_stop_position}")

            except SequencePromptCancelled:

                continue

            except MutationAnalysisError as error:

                style.error(f"{error}")

            except Exception as error:

                style.error(f"Coding analysis failed: {error}")

        # ------------------------------------------------------------------
        # REPORT
        # ------------------------------------------------------------------

        elif choice == "4":

            if last_result is None:

                style.note("Run a sequence comparison first.")

            else:

                print()

                print_mutation_report(
                    last_result,
                    include_alignment=True,
                )

        # ------------------------------------------------------------------
        # SELF TEST
        # ------------------------------------------------------------------

        elif choice == "5":

            tests = (
                run_mutation_self_test()
            )

            style.header('MUTATION SELF-TEST')

            for name, passed in (
                tests.items()
            ):

                print(
                    f"[{'PASS' if passed else 'FAIL'}] "
                    f"{name}"
                )

            print(
                minor_rule()
            )

            print(
                f"{sum(tests.values())}/"
                f"{len(tests)} tests passed."
            )

        else:

            style.invalid_selection()


def run_mutation_variant_ui(
    dna=None,
    header="Sequence",
    sequences=None,
    pipeline=None,
):

    while True:

        style.header('MUTATION & VARIANT ANALYSIS')
        style.field("Sequence", header if dna else "none loaded")
        if dna:
            style.field("Length", f"{len(dna)} bp")

        style.menu([
            ("Analysis", [('1', 'Mutation Analysis'), ('2', 'Variant Analysis')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":

            return run_mutation_analysis_ui(
                dna=dna,
                pipeline=pipeline,
                header=header,
                sequences=sequences,
            )

        if choice == "2":

            try:

                from ui.variant_analysis_menu import (
                    run_variant_analysis_ui
                )

                return run_variant_analysis_ui(
                    pipeline=pipeline,
                    name=(
                        f"Variant Analysis - "
                        f"{header}"
                    ),
                )

            except ImportError as error:

                style.header('VARIANT ANALYSIS IMPORT ERROR')

                print()

                print(
                    f"{type(error).__name__}: "
                    f"{error}"
                )

                print()

                print(
                    "Check that "
                    "analyses/variant_analysis.py "
                    "exists and can be imported."
                )

                style.pause()

            except Exception as error:

                style.header('VARIANT ANALYSIS ERROR')

                print()

                print(
                    f"{type(error).__name__}: "
                    f"{error}"
                )

                style.pause()

        else:

            style.invalid_selection()
