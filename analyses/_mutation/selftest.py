"""The built-in self-test of the mutation module (menu entry and test suite)."""

from __future__ import annotations

from typing import Dict

from analyses._mutation.sequence import (
    iupac_matches,
    normalize_dna,
    reverse_complement,
    reverse_complement_iupac,
    translate_dna,
)
from analyses._mutation.models import MutationAnalysisError
from analyses._mutation.coding import analyze_coding_consequences
from analyses._mutation.comparison import compare_multiple_variants, detect_mutations
from analyses._mutation.context import sequence_quality
from analyses._mutation.impact import (
    crispr_mutation_overlap,
    primer_mutation_overlap,
    restriction_site_impact,
)
from analyses._mutation.reporting import (
    create_advanced_analysis,
    mutation_map,
    result_to_json,
    validate_mutation_result,
)


# ============================================================================
# UI HELPERS
# ============================================================================


# ============================================================================
# MUTATION UI
# ============================================================================


# ============================================================================
# TOP-LEVEL MUTATION + VARIANT UI
# ============================================================================


# ============================================================================
# SELF TEST
# ============================================================================

def run_mutation_self_test() -> Dict[str, bool]:

    tests: Dict[str, bool] = {}

    # ------------------------------------------------------------------
    # SNP
    # ------------------------------------------------------------------

    r = detect_mutations(
        "ATGCCG",
        "ATGTCG",
    )

    tests["snp_detection"] = (
        len(r.mutations) == 1
        and r.mutations[0].mutation_type
        == "SNP"
        and r.mutations[0].position
        == 4
        and r.mutations[0].transition
        is True
    )

    # ------------------------------------------------------------------
    # Multiple substitutions
    # ------------------------------------------------------------------

    r = detect_mutations(
        "ATGCCG",
        "GTGCTG",
    )

    tests["multiple_substitutions"] = (
        len(r.mutations) == 2
    )

    # ------------------------------------------------------------------
    # Insertion
    # ------------------------------------------------------------------

    r = detect_mutations(
        "ATGCCG",
        "ATGGCCG",
    )

    tests["insertion_detection"] = (
        r.statistics.insertions == 1
        and r.statistics.inserted_bases == 1
        and r.statistics.net_length_change == 1
    )

    # ------------------------------------------------------------------
    # Deletion
    # ------------------------------------------------------------------

    r = detect_mutations(
        "ATGGCCG",
        "ATGCCG",
    )

    tests["deletion_detection"] = (
        r.statistics.deletions == 1
        and r.statistics.deleted_bases == 1
        and r.statistics.net_length_change == -1
    )

    # ------------------------------------------------------------------
    # Synonymous
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "GCT",
        "GCC",
    )

    tests["synonymous"] = (
        coding.synonymous_count == 1
        and coding.missense_count == 0
    )

    # ------------------------------------------------------------------
    # Missense
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "GCT",
        "GTT",
    )

    tests["missense"] = (
        coding.missense_count == 1
        and coding.synonymous_count == 0
    )

    # ------------------------------------------------------------------
    # Nonsense
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "CAG",
        "TAG",
    )

    tests["nonsense"] = (
        coding.nonsense_count == 1
    )

    # ------------------------------------------------------------------
    # Start loss
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "ATG",
        "ACG",
    )

    tests["start_loss"] = (
        coding.start_loss_count == 1
    )

    # ------------------------------------------------------------------
    # Stop loss
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "TAA",
        "CAA",
    )

    tests["stop_loss"] = (
        coding.stop_loss_count == 1
    )

    # ------------------------------------------------------------------
    # Frameshift
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "ATGCCCGAA",
        "ATGCCGAA",
    )

    tests["frameshift"] = (
        coding.frameshift is True
    )

    # ------------------------------------------------------------------
    # In-frame indel
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "ATGAAACCC",
        "ATGCCC",
    )

    tests["in_frame_indel"] = (
        coding.frameshift is False
    )

    # ------------------------------------------------------------------
    # IUPAC
    # ------------------------------------------------------------------

    tests["iupac_matching"] = (
        iupac_matches(
            "R",
            "A",
        )
        and not iupac_matches(
            "R",
            "C",
        )
    )

    # ------------------------------------------------------------------
    # IUPAC reverse complement
    # ------------------------------------------------------------------

    tests["iupac_reverse_complement"] = (
        reverse_complement_iupac(
            "ARYN"
        )
        == "NRYT"
    )

    # ------------------------------------------------------------------
    # Translation
    # ------------------------------------------------------------------

    tests["translation"] = (
        translate_dna(
            "ATGGCC"
        )
        == "MA"
    )

    # ------------------------------------------------------------------
    # Reverse complement
    # ------------------------------------------------------------------

    tests["reverse_complement"] = (
        reverse_complement(
            "ATGC"
        )
        == "GCAT"
    )

    # ------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------

    r = detect_mutations(
        "AAGC",
        "GAGG",
    )

    tests["statistics"] = (
        r.statistics.snps == 2
        and r.statistics.transitions == 1
        and r.statistics.transversions == 1
    )

    # ------------------------------------------------------------------
    # Multiple variants
    # ------------------------------------------------------------------

    multi = compare_multiple_variants(
        "ATGCCG",
        [
            "ATGTCG",
            "ATGCCG",
            "ATGGCCG",
        ],
    )

    tests["multiple_variants"] = (
        len(multi) == 3
        and multi[
            "Variant 2"
        ].mutation_count == 0
        and multi[
            "Variant 3"
        ].statistics.insertions == 1
    )

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    tests["serialization"] = (
        '"mutations"'
        in detect_mutations(
            "ATGC",
            "ATGT",
        ).to_json()
    )

    # ------------------------------------------------------------------
    # Advanced analysis
    # ------------------------------------------------------------------

    advanced_result = detect_mutations(
        "ATGCCG",
        "ATGTCG",
    )

    advanced = create_advanced_analysis(
        advanced_result
    )

    tests["advanced_analysis"] = (
        advanced.reference_quality
        is not None
        and advanced.variant_quality
        is not None
        and advanced.spectra
        is not None
        and len(advanced.contexts) == 1
    )

    # ------------------------------------------------------------------
    # JSON advanced
    # ------------------------------------------------------------------

    tests["json_export"] = (
        '"advanced"'
        in result_to_json(
            advanced_result
        )
    )

    # ------------------------------------------------------------------
    # Restriction impact
    # ------------------------------------------------------------------

    restriction = restriction_site_impact(
        "GAATTC",
        "GACTTC",
        {
            "EcoRI": "GAATTC"
        },
    )

    tests["restriction_impact"] = (
        len(restriction) == 1
        and len(
            restriction[0].lost_sites
        ) == 1
    )

    # ------------------------------------------------------------------
    # Primer overlap
    # ------------------------------------------------------------------

    primers = primer_mutation_overlap(
        detect_mutations(
            "ATGCCG",
            "ATGTCG",
        ).mutations,
        [
            ("P1", 1, 5)
        ],
    )

    tests["primer_overlap"] = (
        primers.primer_overlap is True
    )

    # ------------------------------------------------------------------
    # CRISPR overlap
    # ------------------------------------------------------------------

    targets = crispr_mutation_overlap(
        detect_mutations(
            "ATGCCG",
            "ATGTCG",
        ).mutations,
        [
            ("G1", 3, 6)
        ],
    )

    tests["crispr_overlap"] = (
        targets.crispr_target_overlap is True
    )

    # ------------------------------------------------------------------
    # Quality
    # ------------------------------------------------------------------

    quality = sequence_quality(
        "ATGCCG"
    )

    tests["sequence_quality"] = (
        quality.valid
        and quality.length == 6
        and quality.gc_percent
        > 0
    )

    # ------------------------------------------------------------------
    # Invalid sequence
    # ------------------------------------------------------------------

    try:

        normalize_dna(
            "ATGX"
        )

        tests[
            "invalid_dna_rejected"
        ] = False

    except MutationAnalysisError:

        tests[
            "invalid_dna_rejected"
        ] = True

    # ------------------------------------------------------------------
    # Empty sequence
    # ------------------------------------------------------------------

    try:

        normalize_dna(
            ""
        )

        tests[
            "empty_dna_rejected"
        ] = False

    except MutationAnalysisError:

        tests[
            "empty_dna_rejected"
        ] = True

    # ------------------------------------------------------------------
    # Result consistency
    # ------------------------------------------------------------------

    result = detect_mutations(
        "ATGCCG",
        "ATGTCG",
    )

    try:

        tests[
            "result_validation"
        ] = validate_mutation_result(
            result
        )

    except MutationAnalysisError:

        tests[
            "result_validation"
        ] = False

    # ------------------------------------------------------------------
    # Mutation map
    # ------------------------------------------------------------------

    mutation_map_text = mutation_map(
        result.mutations,
        result.reference_length,
    )

    tests["mutation_map"] = (
        "^" in mutation_map_text
    )


    # ------------------------------------------------------------------
    # Professional identity vs. sequence similarity
    # ------------------------------------------------------------------

    deletion_result = detect_mutations(
        "ATGAAACCC",
        "ATGCCC",
    )

    tests["identity_similarity_separation"] = (
        deletion_result.identity_percent == 100.0
        and deletion_result.sequence_similarity_percent
        < deletion_result.identity_percent
        and abs(
            deletion_result.sequence_similarity_percent
            - (100.0 * (1.0 - 3.0 / 9.0))
        ) < 1e-9
    )

    # ------------------------------------------------------------------
    # Ambiguous substitution must not be a complex change
    # ------------------------------------------------------------------

    ambiguous_result = detect_mutations(
        "ATGCCG",
        "ATGRCG",
    )

    tests["ambiguous_substitution_statistics"] = (
        len(ambiguous_result.mutations) == 1
        and ambiguous_result.mutations[0].mutation_type
        == "Substitution"
        and ambiguous_result.statistics.substitutions == 1
        and ambiguous_result.statistics.complex_changes == 0
        and ambiguous_result.statistics.snps == 0
        and ambiguous_result.statistics.unknown_effects == 1
    )

    # ------------------------------------------------------------------
    # TEST_8 regression suite
    # ------------------------------------------------------------------
    # This suite targets the previously observed failures:
    # - nonsense test must contain exactly one SNP
    # - in-frame indels must not create downstream missense cascades
    # - IUPAC R must be accepted
    # - terminal substitutions must remain valid sequence-level events
    test8_ref = (
        "ATGAAACGTGCTGACTTTGGTGACGCTTACCGTAAAGGTCCTTAAACTGCTGACCGTTGACGAA"
        "CTGAAAGCTGGTGCTTACGGTGACCTGATCGCTGAA"
    )

    # Build every test variant from the same reference so that each intended
    # regression contains only the mutation it claims to contain.
    test8_synonymous = test8_ref[:5] + "G" + test8_ref[6:]
    test8_missense = test8_ref[:10] + "T" + test8_ref[11:]

    # Codon 23 is GAA (Glu). GAA -> TAA is a single SNP at position 67.
    test8_nonsense = test8_ref[:66] + "T" + test8_ref[67:]

    # Three independent SNPs: one transition, one transition, one transversion.
    test8_multi = (
        test8_ref[:5] + "G" + test8_ref[6:11]
        + "C" + test8_ref[12:66]
        + "T" + test8_ref[67:]
    )

    # Two transitions: A>G and T>C.
    test8_transitions = (
        test8_ref[:5] + "G" + test8_ref[6:11]
        + "C" + test8_ref[12:]
    )

    # Two transversions: C>A and T>G.
    test8_transversions = (
        test8_ref[:10] + "A" + "G" + test8_ref[12:]
    )

    test8_inframe_insertion = (
        test8_ref[:30] + "GCT" + test8_ref[30:]
    )
    test8_inframe_deletion = (
        test8_ref[:30] + test8_ref[33:]
    )
    test8_frameshift_insertion = (
        test8_ref[:30] + "G" + test8_ref[30:]
    )
    test8_frameshift_deletion = (
        test8_ref[:30] + test8_ref[31:]
    )

    test8_ambiguous = (
        test8_ref[:20] + "R" + test8_ref[21:]
    )

    test8_cases = {
        "TEST_8_SNP_SYNONYMOUS": test8_synonymous,
        "TEST_8_SNP_MISSENSE": test8_missense,
        "TEST_8_SNP_NONSENSE": test8_nonsense,
        "TEST_8_MULTI_SNP": test8_multi,
        "TEST_8_TRANSITIONS": test8_transitions,
        "TEST_8_TRANSVERSIONS": test8_transversions,
        "TEST_8_INFRAME_INSERTION": test8_inframe_insertion,
        "TEST_8_INFRAME_DELETION": test8_inframe_deletion,
        "TEST_8_FRAMESHIFT_INSERTION": test8_frameshift_insertion,
        "TEST_8_FRAMESHIFT_DELETION": test8_frameshift_deletion,
        "TEST_8_AMBIGUOUS": test8_ambiguous,
    }

    try:
        r8 = {
            name: detect_mutations(test8_ref, sequence)
            for name, sequence in test8_cases.items()
        }

        tests["TEST_8_SNP_SYNONYMOUS"] = (
            r8["TEST_8_SNP_SYNONYMOUS"].statistics.snps == 1
            and r8["TEST_8_SNP_SYNONYMOUS"].statistics.synonymous == 1
        )

        tests["TEST_8_SNP_MISSENSE"] = (
            r8["TEST_8_SNP_MISSENSE"].statistics.snps == 1
            and r8["TEST_8_SNP_MISSENSE"].statistics.missense == 1
        )

        tests["TEST_8_SNP_NONSENSE"] = (
            len(r8["TEST_8_SNP_NONSENSE"].mutations) == 1
            and r8["TEST_8_SNP_NONSENSE"].statistics.snps == 1
            and r8["TEST_8_SNP_NONSENSE"].statistics.nonsense == 1
            and r8["TEST_8_SNP_NONSENSE"].statistics.insertions == 0
            and r8["TEST_8_SNP_NONSENSE"].statistics.deletions == 0
            and r8["TEST_8_SNP_NONSENSE"].statistics.frameshifts == 0
        )

        tests["TEST_8_MULTI_SNP"] = (
            r8["TEST_8_MULTI_SNP"].statistics.snps == 3
            and (
                r8["TEST_8_MULTI_SNP"].statistics.transitions
                + r8["TEST_8_MULTI_SNP"].statistics.transversions
            ) == 3
        )

        tests["TEST_8_TRANSITIONS"] = (
            r8["TEST_8_TRANSITIONS"].statistics.snps == 2
            and r8["TEST_8_TRANSITIONS"].statistics.transitions == 2
            and r8["TEST_8_TRANSITIONS"].statistics.transversions == 0
        )

        tests["TEST_8_TRANSVERSIONS"] = (
            r8["TEST_8_TRANSVERSIONS"].statistics.snps == 2
            and r8["TEST_8_TRANSVERSIONS"].statistics.transversions == 2
            and r8["TEST_8_TRANSVERSIONS"].statistics.transitions == 0
        )

        tests["TEST_8_INFRAME_INSERTION"] = (
            r8["TEST_8_INFRAME_INSERTION"].statistics.insertions == 1
            and r8["TEST_8_INFRAME_INSERTION"].statistics.in_frame_indels == 1
            and r8["TEST_8_INFRAME_INSERTION"].statistics.missense == 0
            and r8["TEST_8_INFRAME_INSERTION"].statistics.frameshifts == 0
        )

        tests["TEST_8_INFRAME_DELETION"] = (
            r8["TEST_8_INFRAME_DELETION"].statistics.deletions == 1
            and r8["TEST_8_INFRAME_DELETION"].statistics.in_frame_indels == 1
            and r8["TEST_8_INFRAME_DELETION"].statistics.missense == 0
            and r8["TEST_8_INFRAME_DELETION"].statistics.frameshifts == 0
        )

        tests["TEST_8_FRAMESHIFT_INSERTION"] = (
            r8["TEST_8_FRAMESHIFT_INSERTION"].statistics.insertions == 1
            and r8["TEST_8_FRAMESHIFT_INSERTION"].statistics.frameshifts == 1
        )

        tests["TEST_8_FRAMESHIFT_DELETION"] = (
            r8["TEST_8_FRAMESHIFT_DELETION"].statistics.deletions == 1
            and r8["TEST_8_FRAMESHIFT_DELETION"].statistics.frameshifts == 1
        )

        tests["TEST_8_AMBIGUOUS"] = (
            len(r8["TEST_8_AMBIGUOUS"].mutations) == 1
            and r8["TEST_8_AMBIGUOUS"].mutations[0].effect == "Unknown"
        )

    except Exception:
        for key in test8_cases:
            tests[key] = False

    # ------------------------------------------------------------------
    # Reverse-complement orientation
    # ------------------------------------------------------------------

    reverse_result = detect_mutations(
        test8_ref,
        reverse_complement_iupac(test8_ref),
    )

    tests["TEST_8_REVERSE"] = (
        len(reverse_result.mutations) == 0
        and reverse_result.reference_length
        == reverse_result.variant_length
        and any(
            "reverse complement"
            in note.lower()
            for note in reverse_result.notes
        )
    )

    # ------------------------------------------------------------------
    # Direct coding consequence: reverse-complement orientation
    # ------------------------------------------------------------------

    reverse_coding = analyze_coding_consequences(
        test8_ref,
        reverse_complement_iupac(test8_ref),
    )

    tests["TEST_8_REVERSE_CODING"] = (
        reverse_coding.missense_count == 0
        and reverse_coding.nonsense_count == 0
        and reverse_coding.synonymous_count == 0
        and reverse_coding.start_loss_count == 0
        and reverse_coding.stop_loss_count == 0
        and reverse_coding.frameshift is False
        and reverse_coding.reference_protein
        == reverse_coding.alternative_protein
        and any(
            "orientation-only"
            in note.lower()
            for note in reverse_coding.notes
        )
    )

    # ------------------------------------------------------------------
    # Direct coding consequence: ambiguous IUPAC base
    # ------------------------------------------------------------------

    ambiguous_coding = analyze_coding_consequences(
        test8_ref,
        test8_ambiguous,
    )

    tests["TEST_8_AMBIGUOUS_CODING"] = (
        ambiguous_coding.reference_protein == ""
        and ambiguous_coding.alternative_protein == ""
        and any(
            "ambiguous IUPAC bases"
            in note
            for note in ambiguous_coding.notes
        )
    )

    return tests


def assert_mutation_self_test() -> None:

    results = run_mutation_self_test()

    failed = [
        name
        for name, passed
        in results.items()
        if not passed
    ]

    if failed:

        raise AssertionError(
            "Mutation self-test failed: "
            + ", ".join(failed)
        )
