"""Phase 6 tests: sequence types across modules (audit D14, D13)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tests import _support  # noqa: F401


PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ"
# Every letter is also an IUPAC nucleotide code (audit D14).
IUPAC_LETTER_PEPTIDE = "MKVHDANWYSK" * 4
DNA_SEQ = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"


# ----------------------------------------------------------------------
# Central alphabets and classifier (core/alphabets.py)
# ----------------------------------------------------------------------
class SequenceClassifierTests(unittest.TestCase):
    def classify(self, sequence, declared=None):
        from core import alphabets as A
        return A.classify_sequence(sequence, declared)

    def assertType(self, result, sequence_type, confidence):
        self.assertEqual((result.sequence_type, result.confidence), (sequence_type, confidence),
                         result.describe())

    def test_plain_nucleotides_are_probable_not_confirmed(self):
        self.assertType(self.classify(DNA_SEQ), "DNA", "PROBABLE")
        self.assertType(self.classify(DNA_SEQ.lower()), "DNA", "PROBABLE")
        self.assertType(self.classify("ACGTNNNNACGT"), "DNA", "PROBABLE")
        self.assertType(self.classify(DNA_SEQ.replace("T", "U")), "RNA", "PROBABLE")
        result = self.classify("GGGCCCAAA")  # neither T nor U
        self.assertType(result, "DNA", "PROBABLE")
        self.assertEqual(result.candidates, ("DNA", "RNA", "PROTEIN"))

    def test_protein_only_symbols_confirm_protein(self):
        for sequence in (PROTEIN, "MKXBZUOLLL*", "ACGTE"):
            with self.subTest(sequence=sequence):
                self.assertType(self.classify(sequence), "PROTEIN", "CONFIRMED")

    def test_ambiguity_codes_and_t_u_mixtures_are_ambiguous(self):
        for sequence in (IUPAC_LETTER_PEPTIDE, "ATGCRYATGC", DNA_SEQ[:20] + DNA_SEQ[20:].replace("T", "U")):
            with self.subTest(sequence=sequence):
                result = self.classify(sequence)
                self.assertType(result, None, "AMBIGUOUS")
                self.assertFalse(result.usable)
                self.assertIn("declare", result.reasons[0])
        self.assertEqual(self.classify(IUPAC_LETTER_PEPTIDE).candidates, ("DNA", "RNA", "PROTEIN"))

    def test_invalid_inputs(self):
        for sequence, fragment in (("", "empty"), ("ACGT123", "'1'"), ("ACGT-ACGT", "alignment gaps"),
                                   ("ACGT.ACGT", "alignment gaps")):
            with self.subTest(sequence=sequence):
                result = self.classify(sequence)
                self.assertType(result, None, "INVALID")
                self.assertIn(fragment, result.reasons[0])
        self.assertEqual(self.classify("AC1G2").invalid_symbols, ((3, "1"), (5, "2")))

    def test_declared_type_confirms_or_contradicts(self):
        self.assertType(self.classify(IUPAC_LETTER_PEPTIDE, "protein"), "PROTEIN", "CONFIRMED")
        self.assertType(self.classify("ATGCRYATGC", "DNA"), "DNA", "CONFIRMED")
        self.assertType(self.classify("ACGU", "rna"), "RNA", "CONFIRMED")
        contradicted = self.classify("ACGU", "DNA")
        self.assertType(contradicted, None, "INVALID")
        self.assertIn("declared DNA", contradicted.reasons[0])
        self.assertType(self.classify(PROTEIN, "DNA"), None, "INVALID")
        with self.assertRaises(ValueError):
            self.classify("ACGT", "peptide-ish")

    def test_declared_nucleotide_is_resolved_by_letters(self):
        # GenBank writes mRNA in the DNA alphabet: declared nucleotide + T = DNA letters.
        self.assertType(self.classify("ATGCRY", "NUCLEOTIDE"), "DNA", "CONFIRMED")
        self.assertType(self.classify("AUGCRY", "NUCLEOTIDE"), "RNA", "CONFIRMED")
        self.assertType(self.classify("ATGCU", "NUCLEOTIDE"), None, "INVALID")

    def test_molecule_labels(self):
        from core import alphabets as A
        for label in ("aa", "protein", "Peptide"):
            self.assertEqual(A.declared_type_from_molecule(label), "PROTEIN")
        for label in ("DNA", "mRNA", "ss-RNA", "genomic", "cRNA", "rna"):
            self.assertEqual(A.declared_type_from_molecule(label), "NUCLEOTIDE")
        for label in ("", None, "circular"):
            self.assertIsNone(A.declared_type_from_molecule(label))

    def test_normalization_never_converts_or_drops(self):
        from core import alphabets as A
        self.assertEqual(A.normalize_sequence(" aU g\n-x\t"), "AUG-X")

    def test_require_sequence_type(self):
        from core import alphabets as A
        self.assertEqual(A.require_sequence_type(DNA_SEQ, ("DNA",)).sequence_type, "DNA")
        self.assertEqual(
            A.require_sequence_type(IUPAC_LETTER_PEPTIDE, ("PROTEIN",), declared="PROTEIN").confidence,
            "CONFIRMED")
        for sequence, accepted, fragment in (
            (IUPAC_LETTER_PEPTIDE, ("DNA",), "cannot be determined"),
            (PROTEIN, ("DNA", "RNA"), "needs a DNA or RNA sequence, but this is PROTEIN"),
            ("", ("PROTEIN",), "empty"),
        ):
            with self.subTest(sequence=sequence[:10]):
                with self.assertRaisesRegex(A.SequenceTypeError, fragment) as caught:
                    A.require_sequence_type(sequence, accepted, context="Test")
                self.assertIsNotNone(caught.exception.classification)

    def test_every_nucleotide_code_is_also_an_amino_acid_letter(self):
        """The asymmetry the classifier is built on."""
        from core import alphabets as A
        self.assertLessEqual(A.DNA_IUPAC | A.RNA_IUPAC, A.PROTEIN_IUPAC)
        self.assertEqual(A.PROTEIN_ONLY, frozenset("EFIJLOPQXZ*"))

    @unittest.skipUnless(__import__("importlib").util.find_spec("Bio"), "Biopython (reference extra) not installed")
    def test_alphabets_agree_with_biopython(self):
        """External reference comparison: Bio.Data.IUPACData letter sets."""
        from Bio.Data import IUPACData
        from core import alphabets as A
        self.assertEqual(A.DNA_IUPAC, frozenset(IUPACData.ambiguous_dna_letters))
        self.assertEqual(A.RNA_IUPAC, frozenset(IUPACData.ambiguous_rna_letters))
        self.assertEqual(A.PROTEIN_STANDARD, frozenset(IUPACData.protein_letters))
        self.assertEqual(A.PROTEIN_IUPAC - {"*"}, frozenset(IUPACData.extended_protein_letters))


def _blast_plus_available(*executables: str) -> bool:
    from analyses import blast as B
    manager = B.BlastExecutableManager()
    return all(manager.resolve(name) for name in executables)


# ----------------------------------------------------------------------
# D14-F2: one query-type vocabulary in BLAST
# ----------------------------------------------------------------------
class BlastQueryTypeVocabularyTests(unittest.TestCase):
    def test_program_table_uses_query_vocabulary(self):
        from analyses import blast as B
        for name, info in B.PROGRAMS.items():
            with self.subTest(program=name):
                self.assertIn(info["query_type"], B.QUERY_TYPES)
        self.assertEqual(B.build_query(PROTEIN).sequence_type, "PROTEIN")
        self.assertEqual(B.build_query("ACGTACGT").sequence_type, "DNA")

    def test_programs_offered_for_each_query_type(self):
        from analyses import blast as B
        self.assertEqual(B.programs_for_query_type("PROTEIN"), ["blastp", "tblastn"])
        self.assertEqual(B.programs_for_query_type("DNA"), ["blastn", "blastx", "tblastx"])

    def test_check_accepts_matching_and_rejects_other_type(self):
        from analyses import blast as B
        protein, dna = B.build_query(PROTEIN), B.build_query("ACGTACGTACGT")
        for program in ("blastp", "tblastn"):
            B.check_program_query_type(program, protein)
            with self.assertRaisesRegex(B.BlastAnalysisError, "expects a PROTEIN query"):
                B.check_program_query_type(program, dna)
        for program in ("blastn", "blastx", "tblastx"):
            B.check_program_query_type(program, dna)
            with self.assertRaisesRegex(B.BlastAnalysisError, "expects a DNA query"):
                B.check_program_query_type(program, protein)

    def test_remote_run_rejects_wrong_query_type_before_submitting(self):
        from analyses import blast as B
        client = B.NCBIBlastClient()
        submitted = []
        client.submit = lambda *args, **kwargs: submitted.append(args)  # type: ignore[method-assign]
        with self.assertRaisesRegex(B.BlastAnalysisError, "expects a PROTEIN query"):
            client.run(B.build_query("ACGTACGTACGT"), B.BlastParameters(program="blastp", database="nr"))
        self.assertEqual(submitted, [])

    @unittest.skipUnless(_blast_plus_available("makeblastdb", "blastp"),
                         "BLAST+ (makeblastdb, blastp) not installed")
    def test_local_blastp_runs_protein_query(self):
        """External tool execution: a protein query against a two-record protein database."""
        from analyses import blast as B
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            fasta = Path(tmp) / "proteins.fasta"
            fasta.write_text(f">target1\n{PROTEIN}\n>target2\n{'GSGSWWPPHHCC' * 6}\n")
            database = B.LocalBlastDatabaseManager(root=tmp).create_custom_database(
                fasta, "mini_protein", "protein")
            run = B.LocalBlastEngine().run(
                B.build_query(PROTEIN[5:45], name="Q"),
                B.BlastParameters(program="blastp", database="nr"),
                database,
            )
        self.assertEqual([hit.hit_id for hit in run.hits][:1], ["target1"])


if __name__ == "__main__":
    unittest.main()
