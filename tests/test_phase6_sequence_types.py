"""Phase 6 tests: sequence types across modules (audit D14, D13)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tests import _support  # noqa: F401
from ui import blast_menu


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

    def assertType(self, result, sequence_type, evidence):
        self.assertEqual((result.sequence_type, result.evidence), (sequence_type, evidence),
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
            A.require_sequence_type(IUPAC_LETTER_PEPTIDE, ("PROTEIN",), declared="PROTEIN").evidence,
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
        self.assertEqual((record["alphabet"], record["sequence_classification"]["evidence"]),
                         ("PROTEIN", "CONFIRMED"))

    def test_genbank_mrna_written_with_t_is_a_confirmed_nucleotide(self):
        [record] = self.load("m.gb", _genbank("MRNA1", len(DNA_SEQ), "bp", "mRNA", DNA_SEQ))
        self.assertEqual(record["molecule_type"], "mRNA")
        self.assertEqual(record["declared_type"], "NUCLEOTIDE")
        self.assertEqual((record["alphabet"], record["sequence_classification"]["evidence"]),
                         ("DNA", "CONFIRMED"))

    def test_fasta_moltype_modifier(self):
        records = self.load("x.fasta", (
            f">pep [moltype=protein]\n{IUPAC_LETTER_PEPTIDE}\n"
            f">same_without_modifier\n{IUPAC_LETTER_PEPTIDE}\n"
            f">dna [moltype=genomic DNA]\n{DNA_SEQ}\n"
            f">plain\n{DNA_SEQ}\n"
            f">wrong [moltype=mRNA]\n{PROTEIN}\n"
        ))
        summary = [(r["molecule_type"], r["alphabet"], r["sequence_classification"]["evidence"])
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
        self.assertEqual((sequence_record.alphabet, sequence_record.type_evidence), ("PROTEIN", "CONFIRMED"))
        provenance = sequence_record.provenance()
        self.assertEqual((provenance["sequence_type"], provenance["molecule_type"]), ("PROTEIN", "protein"))

    def test_sequence_record_classifies_and_accepts_a_declaration(self):
        from core.models import SequenceRecord
        self.assertEqual(SequenceRecord("x", "acgt").alphabet, "DNA")
        undecided = SequenceRecord("p", IUPAC_LETTER_PEPTIDE)
        self.assertEqual((undecided.alphabet, undecided.type_evidence), (None, "AMBIGUOUS"))
        declared = undecided.with_declared_type("PROTEIN")
        self.assertEqual((declared.alphabet, declared.type_evidence, declared.declared_type),
                         ("PROTEIN", "CONFIRMED", "PROTEIN"))
        self.assertEqual(SequenceRecord("x", "ACGT", alphabet="DNA").type_evidence, "CONFIRMED")

    def test_report_input_section_names_the_sequence_type(self):
        from core.models import SequenceRecord
        from core.pipeline import SynapsePipeline
        from core.report import build_analysis_report, report_to_markdown
        pipeline = SynapsePipeline(sequence=DNA_SEQ, header="plain")
        pipeline.attach_record(SequenceRecord("plain", DNA_SEQ))
        report = build_analysis_report(pipeline)
        self.assertEqual((report["input"]["sequence_type"], report["input"]["sequence_type_evidence"]),
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
        with mock.patch("builtins.input", return_value="2"), mock.patch("builtins.print"):
            self.assertEqual(blast_menu._query_type_for_menu(IUPAC_LETTER_PEPTIDE, None), "PROTEIN")
        with mock.patch("builtins.input", return_value="0"), mock.patch("builtins.print"):
            self.assertIsNone(blast_menu._query_type_for_menu(IUPAC_LETTER_PEPTIDE, None))
        # A declared or determinable type is not asked for (input() would fail).
        self.assertEqual(blast_menu._query_type_for_menu(IUPAC_LETTER_PEPTIDE, "PROTEIN"), "PROTEIN")
        self.assertEqual(blast_menu._query_type_for_menu(DNA_SEQ, None), "DNA")

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


# ----------------------------------------------------------------------
# D13, D14-F4/F5/F7: one input contract for the menu gate and the API
# ----------------------------------------------------------------------
class InputContractTests(unittest.TestCase):
    def test_both_normalize_dna_functions_share_one_contract(self):
        from core import alphabets, seqtools, validation
        for raw in (" aU g\n", "acgt", "MKV-1"):
            with self.subTest(raw=raw):
                expected = alphabets.normalize_sequence(raw)
                self.assertEqual(validation.normalize_dna(raw), expected)
                self.assertEqual(seqtools.normalize_dna(raw), expected)

    def test_menu_gate_rejects_empty_input(self):
        from core import validation
        self.assertFalse(validation.validate_dna(""))
        self.assertTrue(validation.validate_dna("ACGTN"))
        self.assertEqual(validation.get_dna_errors("AC1U"),
                         [{"character": "1", "position": 3}, {"character": "U", "position": 4}])

    def test_translation_still_reads_rna(self):
        from core.seqtools import get_genetic_code
        self.assertEqual(get_genetic_code(1).translate("AUGGCCUAA"), "MA")

    def test_api_refuses_protein_and_undetermined_input_for_nucleotide_analyses(self):
        from core import api
        from core.alphabets import SequenceTypeError
        calls = (
            lambda s, **k: api.sequence_overview(s, **k),
            lambda s, **k: api.find_orfs(s, **k),
            lambda s, **k: api.analyze_protein_from_cds(s, **k),
            lambda s, **k: api.restriction_digest(s, ["EcoRI"], **k),
        )
        for call in calls:
            with self.subTest(call=call):
                with self.assertRaisesRegex(SequenceTypeError, "but this is PROTEIN"):
                    call(PROTEIN)
                with self.assertRaisesRegex(SequenceTypeError, "cannot be determined"):
                    call(IUPAC_LETTER_PEPTIDE)
                with self.assertRaisesRegex(SequenceTypeError, "empty"):
                    call("")
        # Declaring the type is enough for an undetermined sequence.
        declared = api.sequence_overview(IUPAC_LETTER_PEPTIDE, sequence_type="DNA")
        self.assertEqual(declared.warnings, [])

    def test_api_converts_rna_explicitly_and_records_it(self):
        from core import api
        rna = DNA_SEQ.replace("T", "U")
        for result in (api.sequence_overview(rna), api.find_orfs(rna, min_length=30),
                       api.analyze_protein_from_cds(rna)):
            with self.subTest(module=result.module):
                self.assertIn("RNA input was analyzed in the DNA alphabet (U read as T).", result.warnings)
        self.assertEqual(api.analyze_protein_from_cds(rna).data["protein"],
                         api.analyze_protein_from_cds(DNA_SEQ).data["protein"])
        self.assertEqual(api.sequence_overview(rna).data["gc_content"],
                         api.sequence_overview(DNA_SEQ).data["gc_content"])

    def test_dna_only_analyses_refuse_rna(self):
        from core import api
        from core.alphabets import SequenceTypeError
        rna = DNA_SEQ.replace("T", "U")
        for call in (lambda: api.restriction_digest(rna, ["EcoRI"]), lambda: api.design_crispr_guides(rna)):
            with self.assertRaisesRegex(SequenceTypeError, "needs a DNA sequence, but this is RNA"):
                call()

    def test_api_uses_a_records_declared_type(self):
        from core import api
        from core.models import SequenceRecord
        record = SequenceRecord("pep", IUPAC_LETTER_PEPTIDE, alphabet="PROTEIN")
        self.assertEqual(api.analyze_protein(record).data["length"], len(IUPAC_LETTER_PEPTIDE))
        with self.assertRaisesRegex(ValueError, "but this is PROTEIN"):
            api.find_orfs(record)

    def test_analyze_protein_accepts_proteins_only(self):
        from core import api
        from core.alphabets import SequenceTypeError
        self.assertEqual(api.analyze_protein(PROTEIN).data["length"], len(PROTEIN))
        with self.assertRaisesRegex(SequenceTypeError, "cannot be determined"):
            api.analyze_protein(IUPAC_LETTER_PEPTIDE)
        self.assertEqual(api.analyze_protein(IUPAC_LETTER_PEPTIDE, sequence_type="PROTEIN").module, "protein")
        with self.assertRaisesRegex(SequenceTypeError, "needs a PROTEIN sequence"):
            api.analyze_protein(DNA_SEQ)


# ----------------------------------------------------------------------
# D14-F6: IUPAC protein codes keep their position and meaning
# ----------------------------------------------------------------------
class ProteinResidueCodeTests(unittest.TestCase):
    BASE = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"

    def with_residue(self, residue):
        return self.BASE[:12] + residue + self.BASE[13:]

    def test_every_code_keeps_its_position(self):
        from analyses import protein as P
        for residue in "BZJXUO*?":
            with self.subTest(residue=residue):
                positional = P.get_positional_protein(self.with_residue(residue))
                self.assertEqual(len(positional), len(self.BASE))
                self.assertEqual(positional[12], "?")
                self.assertEqual(positional[13:], self.BASE[13:])
        self.assertEqual(P.get_positional_protein(self.BASE + "*"), self.BASE)

    def test_windows_after_the_code_are_unchanged(self):
        from analyses import protein as P
        clean = {w["start"]: w["score"] for w in P.calculate_hydropathy_profile(self.BASE, 9)}
        coded = {w["start"]: w["score"] for w in P.calculate_hydropathy_profile(self.with_residue("X"), 9)}
        self.assertEqual(len(coded), len(clean) - 9)  # only the 9 windows covering position 13
        for start, score in coded.items():
            self.assertAlmostEqual(score, clean[start])

    def test_quality_categories(self):
        from analyses import protein as P
        for residue, key in (("X", "ambiguous_residues"), ("B", "ambiguous_residues"),
                             ("J", "ambiguous_residues"), ("U", "nonstandard_residues"),
                             ("O", "nonstandard_residues")):
            with self.subTest(residue=residue):
                sequence = self.with_residue(residue)
                quality = P.analyze_protein_quality(sequence)
                self.assertEqual((quality[key], quality["invalid_residues"], quality["valid"]), (1, 0, True))
                self.assertTrue(P.validate_protein_sequence(sequence))
        internal_stop = P.analyze_protein_quality(self.with_residue("*"))
        self.assertEqual((internal_stop["invalid_residues"], internal_stop["valid"]), (1, False))
        terminal = P.analyze_protein_quality(self.BASE + "*")
        self.assertEqual((terminal["length"], terminal["terminal_stop"], terminal["valid"]),
                         (len(self.BASE), True, True))

    def test_molecular_weight_ranges(self):
        from analyses import protein as P
        M = P.AMINO_ACID_RESIDUE_MASSES
        base = P.calculate_protein_molecular_weight(self.BASE)
        without_s = base["min"] - M["S"]
        for residue, (low, high) in (("B", ("D", "N")), ("Z", ("Q", "E")), ("J", ("I", "L"))):
            with self.subTest(residue=residue):
                weight = P.calculate_protein_molecular_weight(self.with_residue(residue))
                expected = sorted((without_s + M[low], without_s + M[high]))
                self.assertAlmostEqual(weight["min"], expected[0])
                self.assertAlmostEqual(weight["max"], expected[1])
                self.assertTrue(weight["ambiguous"])
                self.assertFalse(weight["undetermined"])
        undetermined = P.calculate_protein_molecular_weight(self.with_residue("X"))
        self.assertEqual((undetermined["min"], undetermined["undetermined"]), (None, True))
        self.assertAlmostEqual(P.calculate_protein_molecular_weight(self.BASE + "*")["min"], base["min"])

    @unittest.skipUnless(__import__("importlib").util.find_spec("Bio"), "Biopython (reference extra) not installed")
    def test_sec_and_pyl_masses_agree_with_biopython(self):
        """External reference comparison: Bio.SeqUtils.molecular_weight (average masses)."""
        from Bio.Data.IUPACData import protein_weights
        from Bio.SeqUtils import molecular_weight
        from analyses import protein as P
        for residue in "UO":
            with self.subTest(residue=residue):
                self.assertAlmostEqual(P.AMINO_ACID_RESIDUE_MASSES[residue],
                                       protein_weights[residue] - P.WATER_MASS, places=3)
                sequence = self.with_residue(residue)
                ours = P.calculate_protein_molecular_weight(sequence)["min"]
                # The standard-residue tables (ExPASy vs Biopython) differ by < 0.1 Da here.
                self.assertLess(abs(ours - molecular_weight(sequence, "protein")), 0.1)


# ----------------------------------------------------------------------
# D14-F8, F9, F10: module entry points normalize and check their input
# ----------------------------------------------------------------------
class ModuleInputTests(unittest.TestCase):
    def test_orf_search_is_case_insensitive_and_refuses_other_alphabets(self):
        from analyses import orf as O
        upper = O.find_orfs(DNA_SEQ, min_length=18)
        self.assertGreater(len(upper), 0)
        self.assertEqual(O.find_orfs(DNA_SEQ.lower(), min_length=18), upper)  # was: 0 ORFs
        for sequence in (DNA_SEQ.replace("T", "U"), PROTEIN):
            with self.subTest(sequence=sequence[:10]):
                with self.assertRaisesRegex(ValueError, "ORF search needs a DNA sequence"):
                    O.find_orfs(sequence, min_length=30)

    def test_transcription_checks_its_input(self):
        from analyses.rna import transcribe_rna
        self.assertEqual(transcribe_rna("atg c\n"), "AUGC")
        self.assertEqual(transcribe_rna("AUGC"), "AUGC")
        with self.assertRaisesRegex(ValueError, "found 'E' at"):
            transcribe_rna("MKTAYE")  # was: "MKUAYE"
        with self.assertRaisesRegex(ValueError, "'1' at 4"):
            transcribe_rna("ACG1")

    def test_rna_gc_content_counts_u(self):
        from analyses import sequence as S
        rna = DNA_SEQ.replace("T", "U")
        self.assertAlmostEqual(S.calculate_gc_content(rna), S.calculate_gc_content(DNA_SEQ))
        self.assertEqual(S.count_nucleotides("AUUG")["T"], 2)


class MainMenuSequenceTypeTests(unittest.TestCase):
    """End-to-end: the real main.py with a protein, an undetermined, an RNA and an invalid record."""

    FASTA = (
        f">pep [moltype=protein]\n{IUPAC_LETTER_PEPTIDE}\n"
        ">amb\nATGCRYATGCATGCATGC\n"
        f">rna\n{DNA_SEQ.replace('T', 'U')}\n"
        ">bad\nACGT12ACGT\n"
    )

    def test_each_type_takes_its_route(self):
        import os
        import subprocess
        import sys
        from tests._support import PROJECT_ROOT
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mixed.fasta"
            path.write_text(self.FASTA)
            script = "\n".join([
                str(path),
                "1", "1", "5", "", "0", "S",   # protein record -> protein analysis
                "2", "1", "S",                 # undetermined -> declared DNA -> main menu
                "3", "S",                      # RNA -> main menu in the DNA alphabet
                "4", "",                       # invalid -> validation error
                "1", "0",                      # protein record -> exit
            ]) + "\n"
            env = dict(os.environ, SYNAPSE_LOG_DIR=tmp, PYTHONIOENCODING="utf-8")
            completed = subprocess.run(
                [sys.executable, "main.py"], cwd=PROJECT_ROOT, input=script,
                capture_output=True, text=True, encoding="utf-8", timeout=120, env=env,
            )
        output = completed.stdout + completed.stderr
        self.assertEqual(completed.returncode, 0, output[-2000:])
        self.assertNotIn("Traceback", output)
        self.assertIn(f"pep [moltype=protein]  ({len(IUPAC_LETTER_PEPTIDE)} aa, linear)", output)
        self.assertIn("PROTEIN RECORD", output)
        self.assertIn("COMPLETE PROTEIN REPORT", output.upper())
        self.assertIn("SEQUENCE TYPE", output)
        self.assertIn("RNA sequence, analyzed in the DNA alphabet (U read as T).", output)
        self.assertIn("Invalid character '1' at position 5", output)
        self.assertIn("Thank you for using Synapse.", output)


if __name__ == "__main__":
    unittest.main()
