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
# D14-F3, F11: declared molecule type is read, kept and reported
# ----------------------------------------------------------------------
def _genbank(name, length, unit, molecule, sequence):
    origin = "\n".join(
        f"{i + 1:>9} " + " ".join(sequence[j:j + 10] for j in range(i, min(i + 60, len(sequence)), 10))
        for i in range(0, len(sequence), 60)
    )
    return (f"LOCUS       {name:<16} {length:>11} {unit} {molecule:<6}  linear   SYN 31-OCT-2026\n"
            f"DEFINITION  Synthetic {name}.\nACCESSION   {name}\nORIGIN\n{origin.lower()}\n//\n")


class DeclaredMoleculeTypeTests(unittest.TestCase):
    def load(self, filename, text):
        from fileio.loader import load_sequences
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / filename
            path.write_text(text)
            return load_sequences(str(path))

    def test_genbank_protein_record_is_confirmed_protein(self):
        [record] = self.load("p.gp", _genbank("PEP1", len(IUPAC_LETTER_PEPTIDE), "aa", "", IUPAC_LETTER_PEPTIDE))
        self.assertEqual(record["molecule_type"], "protein")
        self.assertEqual((record["alphabet"], record["sequence_classification"]["confidence"]),
                         ("PROTEIN", "CONFIRMED"))

    def test_genbank_mrna_written_with_t_is_a_confirmed_nucleotide(self):
        [record] = self.load("m.gb", _genbank("MRNA1", len(DNA_SEQ), "bp", "mRNA", DNA_SEQ))
        self.assertEqual(record["molecule_type"], "mRNA")
        self.assertEqual(record["declared_type"], "NUCLEOTIDE")
        self.assertEqual((record["alphabet"], record["sequence_classification"]["confidence"]),
                         ("DNA", "CONFIRMED"))

    def test_fasta_moltype_modifier(self):
        records = self.load("x.fasta", (
            f">pep [moltype=protein]\n{IUPAC_LETTER_PEPTIDE}\n"
            f">same_without_modifier\n{IUPAC_LETTER_PEPTIDE}\n"
            f">dna [moltype=genomic DNA]\n{DNA_SEQ}\n"
            f">plain\n{DNA_SEQ}\n"
            f">wrong [moltype=mRNA]\n{PROTEIN}\n"
        ))
        summary = [(r["molecule_type"], r["alphabet"], r["sequence_classification"]["confidence"])
                   for r in records]
        self.assertEqual(summary, [
            ("protein", "PROTEIN", "CONFIRMED"),
            (None, None, "AMBIGUOUS"),
            ("genomic DNA", "DNA", "CONFIRMED"),
            (None, "DNA", "PROBABLE"),
            ("mRNA", None, "INVALID"),
        ])
        self.assertTrue(any("does not fit the sequence" in note for note in records[4]["notes"]))

    def test_sequence_record_carries_and_exports_the_type(self):
        from core.models import SequenceRecord
        [record] = self.load("x.fasta", f">pep [moltype=protein]\n{IUPAC_LETTER_PEPTIDE}\n")
        sequence_record = SequenceRecord.from_legacy_dict(record, 1)
        self.assertEqual((sequence_record.alphabet, sequence_record.type_confidence), ("PROTEIN", "CONFIRMED"))
        provenance = sequence_record.provenance()
        self.assertEqual((provenance["sequence_type"], provenance["molecule_type"]), ("PROTEIN", "protein"))

    def test_sequence_record_classifies_and_accepts_a_declaration(self):
        from core.models import SequenceRecord
        self.assertEqual(SequenceRecord("x", "acgt").alphabet, "DNA")
        undecided = SequenceRecord("p", IUPAC_LETTER_PEPTIDE)
        self.assertEqual((undecided.alphabet, undecided.type_confidence), (None, "AMBIGUOUS"))
        declared = undecided.with_declared_type("PROTEIN")
        self.assertEqual((declared.alphabet, declared.type_confidence, declared.declared_type),
                         ("PROTEIN", "CONFIRMED", "PROTEIN"))
        self.assertEqual(SequenceRecord("x", "ACGT", alphabet="DNA").type_confidence, "CONFIRMED")

    def test_report_input_section_names_the_sequence_type(self):
        from core.models import SequenceRecord
        from core.pipeline import SynapsePipeline
        from core.report import build_analysis_report, report_to_markdown
        pipeline = SynapsePipeline(sequence=DNA_SEQ, header="plain")
        pipeline.attach_record(SequenceRecord("plain", DNA_SEQ))
        report = build_analysis_report(pipeline)
        self.assertEqual((report["input"]["sequence_type"], report["input"]["sequence_type_confidence"]),
                         ("DNA", "PROBABLE"))
        self.assertEqual(report["input"]["record"]["sequence_type"], "DNA")
        self.assertIn("- sequence_type: `DNA`", report_to_markdown(report))


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



# ----------------------------------------------------------------------
# D14-F1, S20: BLAST takes the query type from core.alphabets
# ----------------------------------------------------------------------
class BlastQueryTypeInferenceTests(unittest.TestCase):
    def test_inference_follows_the_classifier(self):
        from analyses import blast as B
        self.assertEqual(B.infer_sequence_type(DNA_SEQ), "DNA")
        self.assertEqual(B.infer_sequence_type(PROTEIN), "PROTEIN")
        self.assertEqual(B.infer_sequence_type(DNA_SEQ.replace("T", "U")), "RNA")
        self.assertEqual(B.infer_sequence_type(IUPAC_LETTER_PEPTIDE, "PROTEIN"), "PROTEIN")

    def test_undetermined_or_invalid_query_is_refused_not_guessed(self):
        from analyses import blast as B
        # Before: the peptide was run as blastn, RNA and gapped DNA were "invalid protein".
        for sequence, fragment in ((IUPAC_LETTER_PEPTIDE, "specify DNA or PROTEIN"),
                                   ("ATGCATGC" + "AUGCAUGC", "specify DNA or PROTEIN"),
                                   ("ACGT-ACGT", "alignment gaps"),
                                   ("", "empty sequence")):
            with self.subTest(sequence=sequence):
                with self.assertRaisesRegex(B.BlastAnalysisError, fragment):
                    B.build_query(sequence)

    def test_rna_query_is_searched_as_dna_with_a_note(self):
        from analyses import blast as B
        query = B.build_query(DNA_SEQ.replace("T", "U"))
        self.assertEqual((query.sequence_type, query.sequence), ("DNA", DNA_SEQ))
        self.assertEqual(query.notes, ["RNA query searched as DNA: U read as T."])
        self.assertEqual(B.build_query(DNA_SEQ).notes, [])

    def test_protein_alphabet_includes_j_u_o(self):
        from analyses import blast as B
        self.assertEqual(B.build_query("MKTAYJUOAKQRQISF*").sequence_type, "PROTEIN")

    def test_fasta_queries_honor_moltype(self):
        from analyses import blast as B
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "q.fasta"
            path.write_text(f">pep [moltype=protein]\n{IUPAC_LETTER_PEPTIDE}\n>dna\n{DNA_SEQ}\n")
            self.assertEqual([q.sequence_type for q in B.load_fasta_queries(path)], ["PROTEIN", "DNA"])
            path.write_text(f">pep\n{IUPAC_LETTER_PEPTIDE}\n")
            with self.assertRaisesRegex(B.BlastAnalysisError, r"record 1 \('pep'\).*moltype=protein"):
                B.load_fasta_queries(path)

    def test_menu_asks_for_an_undetermined_type(self):
        from unittest import mock
        from analyses import blast as B
        with mock.patch("builtins.input", return_value="2"), mock.patch("builtins.print"):
            self.assertEqual(B._query_type_for_menu(IUPAC_LETTER_PEPTIDE, None), "PROTEIN")
        with mock.patch("builtins.input", return_value="0"), mock.patch("builtins.print"):
            self.assertIsNone(B._query_type_for_menu(IUPAC_LETTER_PEPTIDE, None))
        # A declared or determinable type is not asked for (input() would fail).
        self.assertEqual(B._query_type_for_menu(IUPAC_LETTER_PEPTIDE, "PROTEIN"), "PROTEIN")
        self.assertEqual(B._query_type_for_menu(DNA_SEQ, None), "DNA")

    @unittest.skipUnless(_blast_plus_available("makeblastdb", "blastp"),
                         "BLAST+ (makeblastdb, blastp) not installed")
    def test_declared_iupac_letter_peptide_runs_as_blastp(self):
        """External tool execution: the D14 peptide finds itself once declared a protein."""
        from analyses import blast as B
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            fasta = Path(tmp) / "proteins.fasta"
            fasta.write_text(f">peptide\n{IUPAC_LETTER_PEPTIDE}\n>other\n{PROTEIN}\n")
            database = B.LocalBlastDatabaseManager(root=tmp).create_custom_database(
                fasta, "mini_protein", "protein")
            run = B.LocalBlastEngine().run(
                B.build_query(IUPAC_LETTER_PEPTIDE, sequence_type="PROTEIN"),
                B.BlastParameters(program="blastp", database="nr"),
                database,
            )
        self.assertEqual([hit.hit_id for hit in run.hits][:1], ["peptide"])

if __name__ == "__main__":
    unittest.main()
