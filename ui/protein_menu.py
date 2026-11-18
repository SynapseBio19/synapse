"""Terminal dialogs of analyses/protein.py (moved here in P2-2).

The scientific functions stay in ``analyses.protein``; this module draws the
screens, reads the answers and calls them.
"""

from data.amino_acids import AMINO_ACIDS
from analyses.protein import (  # noqa: F401
    RNA_IUPAC,
    analyze_complete_protein,
    analyze_translation,
    calculate_protein_molecular_weight,
)
from ui import style
from ui.style import major_rule, minor_rule


def _protein_ui_pause():
    style.pause()


def _frame_display(frame):
    if frame is None:
        return "N/A"

    value = str(frame)

    if not value.startswith(("+", "-")):
        value = "+" + value

    return value


def _sequence_applicability(length):
    """
    Qualitative applicability label only.
    It is NOT a statistical confidence probability.
    """
    if length < 15:
        return (
            "LOW",
            "Very short sequence; individual residues can "
            "strongly affect sequence-derived descriptors."
        )

    if length < 30:
        return (
            "MODERATE",
            "Short sequence; sequence-derived descriptors "
            "should be interpreted cautiously."
        )

    return (
        "STANDARD",
        "Sequence length is more suitable for the "
        "displayed sequence-level descriptors."
    )


def _protein_limitations_block(analysis):
    level, reason = _sequence_applicability(
        analysis["length"]
    )

    style.section('ANALYSIS APPLICABILITY')

    style.field('Level', f"{level}")

    style.field('Reason', f"{reason}")

    print()
    print(
        "Applicability level is qualitative; "
        "it is not a statistical confidence probability."
    )

    print(
        "All predictions are sequence-level descriptors "
        "or heuristics."
    )

    print(
        "They do not establish biological function, "
        "experimental stability, 3D structure, "
        "cellular localization or biological activity."
    )


def _protein_input_context(
    analysis,
    source="Selected ORF",
    orf_id=None,
    orf_length=None,
    strand=None,
    reading_frame=None,
    orf_status=None,
):
    style.section('INPUT')

    style.field('Source', f"{source}")

    if orf_id is not None:
        style.field('ORF', f"{orf_id}")

    if orf_length is not None:
        style.field('ORF length', f"{orf_length} bp")

    if strand is not None:
        style.field('Strand', f"{strand}")

    if reading_frame is not None:
        style.field('Reading frame', f"{_frame_display(reading_frame)}")

    if orf_status is not None:
        style.field('ORF status', f"{orf_status}")

    style.field('Protein length', f"{analysis['length']} aa")


def _protein_status_line(analysis):
    quality = analysis["quality"]

    if not quality["valid"]:
        return "INVALID"

    if not quality["fully_resolved"]:
        return "VALID - unresolved residues present"

    return "VALID - fully resolved"


def _protein_quality_block(analysis):
    quality = analysis["quality"]

    style.section('SEQUENCE QUALITY')

    style.field('Standard residues', f"{quality['standard_residues']}")

    if quality.get("nonstandard_residues"):
        print(
            f"Sec/Pyl (U, O):      "
            f"{quality['nonstandard_residues']}"
        )

    style.field('Ambiguous residues', f"{quality['ambiguous_residues']}")

    style.field('Invalid residues', f"{quality['invalid_residues']}")

    style.field('Status', f"{_protein_status_line(analysis)}")


def _protein_composition_block(analysis):
    properties = analysis["properties"]

    style.section('COMPOSITION')

    style.field('Hydrophobic', f"{properties['hydrophobic_percentage']:.2f}%")

    style.field('Polar', f"{properties['polar_percentage']:.2f}%")

    style.field('Positive', f"{properties['positive_percentage']:.2f}%")

    style.field('Negative', f"{properties['negative_percentage']:.2f}%")

    style.field('Aromatic', f"{properties['aromatic_percentage']:.2f}%")


def _protein_key_properties(analysis):
    style.section('KEY PROPERTIES')

    mw = analysis["molecular_weight"]

    if mw["ambiguous"]:
        style.field('Molecular weight', f"{mw['min']:,.0f} - "
            f"{mw['max']:,.0f} Da")
    else:
        style.field('Molecular weight', f"{mw['min']:,.0f} Da")

    if analysis["pI"] is not None:
        style.field('Theoretical pI', f"{analysis['pI']:.2f}")

    if analysis["gravy"] is not None:
        style.field('GRAVY', f"{analysis['gravy']:+.3f}")

    if analysis["instability_index"] is not None:
        style.field('Instability index', f"{analysis['instability_index']:.2f}")

    style.field('Aliphatic index', f"{analysis['aliphatic_index']:.2f}")

    style.field('Aromaticity', f"{analysis['aromaticity']:.3f}")


def _translation_quality_block(translation):
    style.section('TRANSLATION QUALITY')

    style.field('Start codon', f"{translation['start_codon'] or 'Not detected'}")

    style.field('Start detected', f"{'Yes' if translation['start_codon'] else 'No'}")

    style.field('Stop codon', f"{translation['stop_codon'] if translation['stop_found'] else 'Not detected'}")

    style.field('Complete endpoint', f"{'Yes' if translation['stop_found'] else 'No'}")

    print(
        f"Trailing nucleotides:"
        f" {translation['trailing_nucleotides']}"
    )

    style.field('Translation status', f"{'Complete' if translation['stop_found'] else 'Partial'}")


def _safe_extinction_values(analysis):
    extinction = analysis["extinction_coefficient"]

    reduced = extinction.get(
        "reduced",
        extinction.get(
            "reduced_coefficient",
            extinction.get(
                "reduced_cysteine",
                0
            )
        )
    )

    oxidized = extinction.get(
        "oxidized",
        extinction.get(
            "oxidized_coefficient",
            extinction.get(
                "oxidized_cysteine",
                0
            )
        )
    )

    return (
        extinction,
        reduced,
        oxidized
    )


def _hydropathy_interpretation(gravy):
    if gravy is None:
        return "Unavailable"

    if gravy >= 1.0:
        return "Strong hydrophobic tendency"

    if gravy >= 0.3:
        return "Moderate hydrophobic tendency"

    if gravy <= -1.0:
        return "Strong hydrophilic tendency"

    if gravy < -0.3:
        return "Hydrophilic tendency"

    return "Mixed hydrophobic / hydrophilic character"


def _charge_interpretation(analysis):
    charges = analysis.get(
        "charges",
        {}
    )

    charge = charges.get(
        "7.4"
    )

    if charge is None:
        charge = charges.get(
            "7.0"
        )

    if charge is None:
        return "Unavailable"

    if charge > 1:
        return "Net positive tendency near neutral pH"

    if charge < -1:
        return "Net negative tendency near neutral pH"

    return "Approximately neutral near neutral pH"


def _stability_interpretation(analysis):
    index = analysis.get(
        "instability_index"
    )

    if index is None:
        return "Unavailable - fewer than 2 resolved residues"

    if index < 40:
        return (
            "Sequence-based stable tendency "
            "(II < 40)"
        )

    return (
        "Sequence-based unstable tendency "
        "(II >= 40)"
    )


def _print_protein_overview(
    analysis,
    source="Selected ORF",
    orf_id=None,
    orf_length=None,
    strand=None,
    reading_frame=None,
    orf_status=None,
):
    style.header(
        "PROTEIN OVERVIEW"
    )

    _protein_input_context(
        analysis,
        source=source,
        orf_id=orf_id,
        orf_length=orf_length,
        strand=strand,
        reading_frame=reading_frame,
        orf_status=orf_status,
    )

    style.section('PROTEIN')

    style.field('Length', f"{analysis['length']} aa")

    style.field('Sequence', f"{analysis['protein']}")

    translation = analysis["translation"]

    if translation is not None:
        style.field('Translation', f"{'Complete' if translation['stop_found'] else 'Partial'}")

    print()
    _protein_quality_block(
        analysis
    )

    print()
    _protein_key_properties(
        analysis
    )

    print()
    _protein_composition_block(
        analysis
    )

    style.section('ANALYSIS QUALITY')

    level, reason = _sequence_applicability(
        analysis["length"]
    )

    style.field('Applicability level', f"{level}")

    style.field('Reason', f"{reason}")

    if analysis["length"] < 30:
        print()
        print(
            "This is a qualitative applicability label, "
            "not a statistical confidence probability."
        )

    print()
    print(
        "This overview contains sequence-level "
        "descriptors only."
    )


def _print_translation_and_sequence(analysis):
    style.header(
        "TRANSLATION & SEQUENCE"
    )

    translation = analysis["translation"]

    if translation is None:
        print(
            "Translation context: Not supplied."
        )

        print(
            "Protein sequence was supplied directly."
        )

    else:
        style.section('TRANSLATION')

        style.field('Reading frame', f"{_frame_display(translation['reading_frame'])}")

        style.field('First codon', f"{translation['start_codon'] or 'N/A'}")

        style.field('AUG at start', f"{'Yes' if translation['starts_with_aug'] else 'No'}")

        style.field('Codons processed', f"{translation['codons_processed']}")

        _translation_quality_block(
            translation
        )

        style.section('CODON TRANSLATION')

        for record in translation[
            "codon_records"
        ]:
            aa = record[
                "amino_acid"
            ]

            if record[
                "status"
            ] == "stop":
                result = "STOP"

            elif aa is not None:
                result = aa

            elif record[
                "possibilities"
            ]:
                result = "/".join(
                    record[
                        "possibilities"
                    ]
                )

            else:
                result = "?"

            print(
                f"{record['number']:>3}: "
                f"{record['codon']} -> "
                f"{result}"
            )

    style.section('PROTEIN SEQUENCE')

    sequence = analysis[
        "protein"
    ]

    width = 60

    for start in range(
        0,
        len(sequence),
        width
    ):
        end = min(
            start + width,
            len(sequence)
        )

        print(
            f"{start + 1:>5}-{end:<5} "
            f"{sequence[start:end]}"
        )

    print()
    style.field('Length', f"{len(sequence)} aa")

    print()
    _protein_limitations_block(
        analysis
    )


def _print_amino_acid_composition(analysis):
    style.header(
        "AMINO ACID COMPOSITION"
    )

    composition = analysis[
        "composition"
    ]

    percentages = analysis[
        "percentages"
    ]

    counts = composition[
        "counts"
    ]

    widths = (3, 4, 16, 11, 6)
    print(style.columns(("AA", "3L", "Name", "Type", "Count", "Percent"), widths))

    print(minor_rule())

    ordered = sorted(
        counts.items(),
        key=lambda item: (
            -item[1],
            item[0]
        )
    )

    for code, count in ordered:
        amino_acid = AMINO_ACIDS[
            code
        ]

        print(style.columns((
            code, amino_acid["three_letter"], amino_acid["name"],
            amino_acid["type"], count, f"{percentages[code]:.2f}%",
        ), widths)
        )

    print(minor_rule())

    style.field('Resolved residues', f"{composition['resolved_length']}")

    style.field('Unresolved residues', f"{composition['ambiguous_count']}")

    style.section('COMPOSITION GROUPS')

    _protein_composition_block(
        analysis
    )

    if counts:
        residue, count = max(
            counts.items(),
            key=lambda item: item[1]
        )

        name = AMINO_ACIDS[
            residue
        ][
            "name"
        ]

        percent = percentages[
            residue
        ]

        print()
        style.field('Most abundant residue', f"{residue} ({name})")

        style.field('Occurrences', f"{count} / "
            f"{composition['resolved_length']}")

        style.field('Relative abundance', f"{percent:.2f}%")

    print()
    _protein_limitations_block(
        analysis
    )


def _print_physicochemical_structural(analysis):
    style.header(
        "PHYSICOCHEMICAL & STRUCTURAL PROPERTIES"
    )

    print()
    _protein_limitations_block(
        analysis
    )

    style.section('MOLECULAR PROPERTIES')

    _protein_key_properties(
        analysis
    )

    print()
    style.field('Hydropathy', f"{_hydropathy_interpretation(analysis['gravy'])}")

    style.field('Charge', f"{_charge_interpretation(analysis)}")

    style.field('Stability', f"{_stability_interpretation(analysis)}")

    style.section('CHARGE PROFILE')

    for ph, charge in analysis[
        "charges"
    ].items():
        print(
            f"Net charge pH {ph:>4}: "
            f"{charge:+.2f}"
        )

    style.section('EXTINCTION COEFFICIENT')

    extinction, reduced, oxidized = (
        _safe_extinction_values(
            analysis
        )
    )

    style.field('Wavelength', "280 nm")

    style.field('Tryptophan', f"{extinction['tryptophan']}")

    style.field('Tyrosine', f"{extinction['tyrosine']}")

    style.field('Cysteine', f"{extinction['cysteine']}")

    style.field('All Cys reduced', f"{reduced} M^-1 cm^-1")

    style.field('Max. disulfide form', f"{oxidized} M^-1 cm^-1")

    style.field('Disulfide pairs', f"{extinction['disulfide_pairs']}")

    if extinction["unpaired_cysteine"]:
        style.field('Unpaired cysteine', "1")

    style.section('SECONDARY STRUCTURE SCREENING')

    secondary = analysis[
        "secondary_structure"
    ]

    style.field('Helix propensity', f"{secondary['helix']:.3f}")

    style.field('Turn propensity', f"{secondary['turn']:.3f}")

    style.field('Sheet propensity', f"{secondary['sheet']:.3f}")

    style.field('Interpretation', "Chou-Fasman sequence propensity")

    style.field('Limitation', "Not a 3D structure prediction")

    style.section('HYDROPHOBIC REGION SCREENING')

    regions = analysis[
        "hydrophobic_regions"
    ]

    if regions:
        for index, region in enumerate(
            regions,
            start=1
        ):
            print(
                f"Region {index}:        "
                f"{region['start']}-"
                f"{region['end']} | "
                f"Hydropathy "
                f"{region['hydropathy']:+.2f}"
            )
    else:
        print(
            "No strongly hydrophobic windows detected."
        )

    style.section('TRANSMEMBRANE SCREENING')

    tm_domains = analysis[
        "transmembrane_domains"
    ]

    if tm_domains:
        for index, domain in enumerate(
            tm_domains,
            start=1
        ):
            print(
                f"TM {index}:             "
                f"{domain['start']}-"
                f"{domain['end']} | "
                f"Hydropathy "
                f"{domain['hydropathy']:+.2f}"
            )
    else:
        style.field('Prediction', "No strong TM-like region detected")

    style.field('Method', "Sequence-only heuristic")

    style.section('SIGNAL PEPTIDE SCREENING')

    signal = analysis[
        "signal_peptide"
    ]

    if signal[
        "predicted"
    ]:
        style.field('Prediction', "Possible N-terminal signal peptide")

        style.field('Heuristic evidence', f"{signal['evidence']}")
    else:
        style.field('Prediction', "No convincing pattern")

    style.field('Method', "Sequence-only heuristic")

    style.section('FLEXIBILITY SCREENING')

    flexibility = analysis[
        "flexibility"
    ]

    if flexibility:
        scores = [
            entry["score"]
            for entry in flexibility
        ]

        style.field('Average', f"{sum(scores) / len(scores):.3f}")

        style.field('Minimum', f"{min(scores):.3f}")

        style.field('Maximum', f"{max(scores):.3f}")
    else:
        print(
            "Sequence too short for "
            "window-based analysis."
        )

    print()
    print(
        "Structural results are sequence-level "
        "screening outputs, not experimental "
        "structure determinations."
    )


def _print_complete_protein_report(
    analysis,
    source="Selected ORF",
    orf_id=None,
    orf_length=None,
    strand=None,
    reading_frame=None,
    orf_status=None,
):
    style.header(
        "COMPLETE PROTEIN REPORT"
    )

    style.section('1. INPUT')

    style.field('Source', f"{source}")

    if orf_id is not None:
        style.field('ORF', f"{orf_id}")

    if orf_length is not None:
        style.field('ORF length', f"{orf_length} bp")

    if strand is not None:
        style.field('Strand', f"{strand}")

    if reading_frame is not None:
        style.field('Reading frame', f"{_frame_display(reading_frame)}")

    if orf_status is not None:
        style.field('ORF status', f"{orf_status}")

    style.field('Protein length', f"{analysis['length']} aa")

    style.section('2. TRANSLATION')

    translation = analysis[
        "translation"
    ]

    if translation is None:
        print(
            "Translation context: Not supplied"
        )
    else:
        style.field('Reading frame', f"{_frame_display(translation['reading_frame'])}")

        style.field('Start codon', f"{translation['start_codon'] or 'Not detected'}")

        style.field('Stop codon', f"{translation['stop_codon'] if translation['stop_found'] else 'Not detected'}")

        style.field('Complete endpoint', f"{'Yes' if translation['stop_found'] else 'No'}")

        style.field('Translation status', f"{'Complete' if translation['stop_found'] else 'Partial'}")

        print(
            f"Trailing nucleotides:"
            f" {translation['trailing_nucleotides']}"
        )

    style.field('Protein sequence', f"{analysis['protein']}")

    style.section('3. SEQUENCE QUALITY')

    _protein_quality_block(
        analysis
    )

    style.section('4. COMPOSITION')

    _protein_composition_block(
        analysis
    )

    counts = analysis[
        "composition"
    ][
        "counts"
    ]

    if counts:
        residue, count = max(
            counts.items(),
            key=lambda item: item[1]
        )

        style.field('Most abundant residue', f"{residue} "
            f"({AMINO_ACIDS[residue]['name']})")

        style.field('Occurrences', f"{count}")

    style.section('5. PHYSICOCHEMISTRY')

    _protein_key_properties(
        analysis
    )

    style.field('Hydropathy', f"{_hydropathy_interpretation(analysis['gravy'])}")

    style.field('Charge', f"{_charge_interpretation(analysis)}")

    style.field('Stability', f"{_stability_interpretation(analysis)}")

    style.section('6. STRUCTURAL SCREENING')

    secondary = analysis[
        "secondary_structure"
    ]

    style.field('Helix propensity', f"{secondary['helix']:.3f}")

    style.field('Turn propensity', f"{secondary['turn']:.3f}")

    style.field('Sheet propensity', f"{secondary['sheet']:.3f}")

    style.field('TM-like regions', f"{len(analysis['transmembrane_domains'])}")

    signal = analysis[
        "signal_peptide"
    ]

    style.field('Signal peptide', f"{'Possible' if signal['predicted'] else 'Not detected'}")

    style.field('Structural method', "Sequence-level heuristics")

    style.section('7. EXTINCTION COEFFICIENT')

    extinction, reduced, oxidized = (
        _safe_extinction_values(
            analysis
        )
    )

    style.field('Wavelength', "280 nm")

    style.field('Tryptophan', f"{extinction['tryptophan']}")

    style.field('Tyrosine', f"{extinction['tyrosine']}")

    style.field('Cysteine', f"{extinction['cysteine']}")

    style.field('All Cys reduced', f"{reduced} M^-1 cm^-1")

    style.field('Max. disulfide form', f"{oxidized} M^-1 cm^-1")

    style.field('Disulfide pairs', f"{extinction['disulfide_pairs']}")

    if extinction["unpaired_cysteine"]:
        style.field('Unpaired cysteine', "1")

    style.section('8. INTERPRETATION')

    level, reason = _sequence_applicability(
        analysis["length"]
    )

    style.field('Applicability level', f"{level}")

    style.field('Reason', f"{reason}")

    print()
    print(
        "All predictions are sequence-level "
        "descriptors or heuristics."
    )

    print(
        "They do not establish biological function, "
        "experimental stability, 3D structure, "
        "cellular localization or biological activity."
    )

    print()
    print(major_rule())


def _protein_ui_menu():

    style.header('PROTEIN ANALYSIS')

    style.menu([
        ("Overview", [('1', 'Protein Overview'), ('2', 'Translation & Sequence')]),
        ("Analysis", [('3', 'Amino Acid Composition'),
                      ('4', 'Physicochemical & Structural Properties')]),
        ("Report", [('5', 'Complete Protein Report')]),
    ])


def protein_analysis_menu(
    protein,
    rna=None,
    source="Selected ORF",
    orf_id=None,
    orf_length=None,
    strand=None,
    reading_frame=None,
    orf_status=None,
    initiator_start=False,
    genetic_code=1,
):
    """
    Main Protein Analysis UI.

    If RNA is supplied, it is used for translation.
    There is no complete-DNA fallback.
    """

    if rna is not None:

        translation = analyze_translation(
            rna,
            initiator_start=initiator_start,
            genetic_code=genetic_code,
        )

        # Preserve the biological frame from the selected ORF.
        if reading_frame is not None:
            translation[
                "reading_frame"
            ] = reading_frame

        active_protein = translation[
            "protein"
        ]

    else:

        translation = None
        active_protein = protein

    if not isinstance(
        active_protein,
        str
    ):
        active_protein = ""

    active_protein = (
        active_protein
        .upper()
    )

    if not active_protein:

        style.header(
            "PROTEIN ANALYSIS ERROR"
        )

        style.note("No resolvable protein sequence is available.")

        _protein_ui_pause()

        return

    analysis = analyze_complete_protein(
        active_protein,
        rna=None,
    )

    # Recalculate molecular weight with the actual translation
    # records so ambiguous codons can contribute a mass range.
    if translation is not None:
        analysis[
            "molecular_weight"
        ] = calculate_protein_molecular_weight(
            active_protein,
            translation
        )

    analysis[
        "translation"
    ] = translation

    while True:

        _protein_ui_menu()

        selection = style.ask_choice()

        if selection == "0":
            return

        elif selection == "1":

            _print_protein_overview(
                analysis,
                source=source,
                orf_id=orf_id,
                orf_length=orf_length,
                strand=strand,
                reading_frame=reading_frame,
                orf_status=orf_status,
            )

            _protein_ui_pause()

        elif selection == "2":

            _print_translation_and_sequence(
                analysis
            )

            _protein_ui_pause()

        elif selection == "3":

            _print_amino_acid_composition(
                analysis
            )

            _protein_ui_pause()

        elif selection == "4":

            _print_physicochemical_structural(
                analysis
            )

            _protein_ui_pause()

        elif selection == "5":

            _print_complete_protein_report(
                analysis,
                source=source,
                orf_id=orf_id,
                orf_length=orf_length,
                strand=strand,
                reading_frame=reading_frame,
                orf_status=orf_status,
            )

            _protein_ui_pause()

        else:

            style.invalid_selection()


def run_protein_analysis_from_record(
    protein,
    name="Protein record",
    pipeline=None,
):
    """
    Start Protein Analysis from a protein sequence (a protein FASTA or
    GenPept record; audit D14).

    The sequence must be in the IUPAC protein alphabet. When a pipeline is
    given, the complete protein analysis is stored before the interactive
    menu opens.
    """

    from core.alphabets import PROTEIN as _PROTEIN, require_sequence_type

    require_sequence_type(
        protein,
        (_PROTEIN,),
        declared=_PROTEIN,
        context="Protein analysis",
    )

    clean_protein = "".join(
        protein.upper().split()
    )

    if pipeline is not None:
        context = getattr(pipeline, "context", None)
        analysis = analyze_complete_protein(clean_protein)
        if context is not None:
            context.protein_source = "Protein record"
            context.store_protein_analysis(analysis, protein_sequence=clean_protein)

    protein_analysis_menu(
        protein=clean_protein,
        source="Protein record",
        orf_id=name,
    )


def run_protein_analysis_from_orf(
    orf_sequence,
    orf_id="Selected ORF",
    strand=None,
    reading_frame=None,
    orf_status=None,
    pipeline=None,
    genetic_code=1,
):
    """
    Start Protein Analysis from exactly one selected ORF.

    When a pipeline is given, the complete protein analysis (computed with the
    initiator-Met rule) is stored before the interactive menu opens.

    The complete DNA sequence is never used as a fallback.
    """

    if not isinstance(
        orf_sequence,
        str
    ):

        print(
            "Protein analysis error: "
            "selected ORF sequence must be a string."
        )

        return

    clean_orf = "".join(
        orf_sequence
        .upper()
        .split()
    )

    allowed = set(
        RNA_IUPAC
    ) | {"T"}

    if (
        not clean_orf
        or any(
            base not in allowed
            for base in clean_orf
        )
    ):

        print(
            "Protein analysis error: "
            "selected ORF contains invalid "
            "nucleotide symbols."
        )

        return

    rna = clean_orf.replace(
        "T",
        "U"
    )

    # A Synapse ORF begins at its start codon, so the first codon is decoded
    # as initiator Met (also for GTG/TTG starts) - except a 5'-partial ORF,
    # whose start lies upstream of the sequence.
    has_start = orf_status not in ("Partial - Missing Start", "Partial - Missing Start and Stop")
    translation = analyze_translation(
        rna,
        initiator_start=has_start,
        genetic_code=genetic_code,
    )

    if reading_frame is not None:

        translation[
            "reading_frame"
        ] = reading_frame

    protein = translation[
        "protein"
    ]

    if not protein:

        print(
            "Protein analysis error: "
            "selected ORF produced no "
            "resolvable protein sequence."
        )

        return

    if pipeline is not None:
        context = getattr(pipeline, "context", None)
        try:
            analysis = analyze_complete_protein("", rna=rna, initiator_start=has_start, genetic_code=genetic_code)
            if context is not None:
                context.protein_source = "Selected ORF"
                context.store_protein_analysis(analysis, protein_sequence=analysis.get("protein", ""))
        except Exception:
            import logging
            logging.getLogger("synapse.protein").exception(
                "Could not store protein result in the pipeline."
            )

    protein_analysis_menu(
        protein=protein,
        rna=rna,
        source="Selected ORF",
        orf_id=orf_id,
        orf_length=len(clean_orf),
        strand=strand,
        reading_frame=(
            reading_frame
            if reading_frame is not None
            else "+1"
        ),
        orf_status=orf_status,
        initiator_start=has_start,
        genetic_code=genetic_code,
    )
