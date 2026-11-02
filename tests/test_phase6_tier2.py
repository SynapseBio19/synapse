"""Phase 6 Tier 2 tests: functional defects (docs/PHASE6_INITIAL_AUDIT.md)."""

from __future__ import annotations

import pathlib
import unittest

from tests import _support  # noqa: F401
from tests._support import PROJECT_ROOT

LOCAL_TEST_DB = PROJECT_ROOT / "data" / "local_test" / "synapse_test"


def _blast_plus_available(*executables: str) -> bool:
    from analyses import blast as B
    manager = B.BlastExecutableManager()
    return all(manager.resolve(name) for name in executables)


# ----------------------------------------------------------------------
# D6, S7: programmatic local BLAST with a database path prefix
# ----------------------------------------------------------------------
class LocalBlastDatabasePathTests(unittest.TestCase):
    def _query(self):
        lines = (PROJECT_ROOT / "data" / "local_test" / "query.fasta").read_text().splitlines()
        return "".join(line.strip() for line in lines if not line.startswith(">"))

    @unittest.skipUnless(_blast_plus_available("blastn", "blastdbcmd"),
                         "BLAST+ (blastn, blastdbcmd) not installed")
    def test_prefix_path_runs_and_records_provenance(self):
        """External tool execution: the bundled database addressed by its prefix."""
        from analyses import blast as B
        run = B.run_blast_analysis(self._query(), engine="local", program="blastn",
                                   local_database_path=str(LOCAL_TEST_DB))
        self.assertEqual([(hit.hit_id, round(hit.best_coverage_percent, 2)) for hit in run.hits],
                         [("Test_Gene_001", 100.0), ("Test_Gene_002", 98.04)])  # audit section 6
        details = run.metadata.database_details
        self.assertEqual(details["title"], "Synapse Local Test Database")
        self.assertEqual(details["sequence_count"], 3)
        self.assertIsNotNone(details["database_date"])

    def test_missing_prefix_is_reported(self):
        from analyses import blast as B
        with self.assertRaisesRegex(B.BlastAnalysisError, "No BLAST database found at prefix"):
            B.run_blast_analysis(self._query(), engine="local", program="blastn",
                                 local_database_path=str(LOCAL_TEST_DB.parent / "missing_db"))

    @unittest.skipUnless(_blast_plus_available("blastdbcmd"), "BLAST+ (blastdbcmd) not installed")
    def test_program_database_type_mismatch_is_refused(self):
        from analyses import blast as B
        with self.assertRaisesRegex(B.BlastAnalysisError, "requires a protein database"):
            B.run_blast_analysis("MKTAYIAKQRQISFVKSHFSRQ", engine="local", program="blastp",
                                 local_database_path=str(LOCAL_TEST_DB))


# ----------------------------------------------------------------------
# D9: the guide's own site is not an off-target
# ----------------------------------------------------------------------
def _random_dna(length, seed):
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


class CrisprOnTargetTests(unittest.TestCase):
    SEQUENCE = _random_dna(300, 7)

    def test_own_site_is_listed_but_not_penalized(self):
        from analyses import crispr as C
        result = C.design_guides(self.SEQUENCE, local_reference_sequences={"self": self.SEQUENCE})
        self.assertTrue(result.candidates)
        self.assertEqual({c.strand for c in result.candidates}, {"+", "-"})
        for candidate in result.candidates:
            with self.subTest(candidate=candidate.candidate_id):
                profile = candidate.specificity_profile
                self.assertEqual((profile["on_target_matches"], profile["exact_matches"]), (1, 0))
                self.assertEqual(candidate.off_target_hits, profile["total_hits"])
                own = [h for h in candidate.local_off_target_hits if h["on_target"]]
                self.assertEqual(len(own), 1)
                self.assertEqual(own[0]["classification"], "on-target")
                self.assertEqual((own[0]["guide_start"], own[0]["strand"]),
                                 (candidate.guide_start, candidate.strand))
        best = result.candidates[0]
        self.assertEqual(best.rank_components["Specificity"], 100.0)  # was 60.0 (100 - 40)
        self.assertNotIn("INFO: local-sequence matches found", best.flags)

    def test_an_identical_copy_elsewhere_is_still_an_off_target(self):
        from analyses import crispr as C
        result = C.design_guides(self.SEQUENCE, local_reference_sequences={
            "self": self.SEQUENCE, "copy": self.SEQUENCE})
        best = result.candidates[0]
        self.assertEqual((best.specificity_profile["on_target_matches"],
                          best.specificity_profile["exact_matches"]), (1, 1))
        self.assertEqual(best.rank_components["Specificity"], 60.0)

    def test_near_match_in_another_sequence_is_counted(self):
        from analyses import crispr as C
        result = C.design_guides(self.SEQUENCE)
        candidate = result.candidates[0]
        mutated = candidate.guide[:2] + ("A" if candidate.guide[2] != "A" else "C") + candidate.guide[3:]
        other = _random_dna(40, 1) + mutated + candidate.pam + _random_dna(40, 2)
        C.annotate_local_off_targets([candidate], {"self": self.SEQUENCE, "other": other})
        profile = candidate.specificity_profile
        self.assertEqual(profile["on_target_matches"], 1)
        self.assertGreaterEqual(profile["one_mismatch_hits"], 1)


# ----------------------------------------------------------------------
# D8: modeled cut site, from cited conventions
# ----------------------------------------------------------------------
class CrisprCutSiteTests(unittest.TestCase):
    SEQUENCE = _random_dna(400, 3)

    def _check(self, system):
        from analyses import crispr as C
        from core.seqtools import reverse_complement
        definition = C.get_crispr_system(system)
        result = C.design_guides(self.SEQUENCE, system=system)
        strands = set()
        for candidate in result.candidates:
            strands.add(candidate.strand)
            for site, k in ((candidate.cut_site, definition.cut_after),
                            (candidate.cut_site_complementary, definition.cut_after_complementary)):
                if candidate.strand == "+":
                    before, after = self.SEQUENCE[:site], self.SEQUENCE[site:]
                else:
                    before, after = reverse_complement(self.SEQUENCE[site:]), reverse_complement(self.SEQUENCE[:site])
                self.assertTrue(before.endswith(candidate.guide[:k]), candidate.candidate_id)
                self.assertTrue(after.startswith(candidate.guide[k:]), candidate.candidate_id)
        self.assertEqual(strands, {"+", "-"})
        return result

    def test_cas9_cuts_blunt_three_nt_from_the_pam(self):
        for system in ("SpCas9", "SpCas9-NG", "SaCas9"):
            with self.subTest(system=system):
                for candidate in self._check(system).candidates:
                    self.assertEqual(candidate.cut_type, "blunt")
                    self.assertEqual(candidate.cut_site, candidate.cut_site_complementary)
                    if candidate.strand == "+":
                        self.assertEqual(candidate.pam_start - (candidate.cut_site + 1), 3)
                    else:
                        self.assertEqual(candidate.cut_site - candidate.pam_end, 3)

    def test_cas12a_cuts_staggered_with_a_five_nt_overhang(self):
        for candidate in self._check("Cas12a").candidates:
            self.assertEqual(candidate.cut_type, "staggered, 5-nt 5' overhang")
            self.assertEqual(abs(candidate.cut_site_complementary - candidate.cut_site), 5)

    def test_cut_is_reported_and_exported(self):
        import json
        from analyses import crispr as C
        result = C.design_guides(self.SEQUENCE)
        candidate = result.candidates[0]
        self.assertIn(f"between {candidate.cut_site} and {candidate.cut_site + 1}",
                      C.format_candidate(candidate, 1))
        exported = json.loads(C.result_to_json(result))["candidates"][0]
        self.assertEqual((exported["cut_site"], exported["cut_type"]), (candidate.cut_site, "blunt"))
        self.assertIn("PMID 22745249", C.format_system_details("SpCas9"))


# ----------------------------------------------------------------------
# D10, D11: CRISPR configuration fields act; exports stay in the project
# ----------------------------------------------------------------------
class CrisprConfigurationTests(unittest.TestCase):
    SEQUENCE = _random_dna(300, 11) + "A" * 9 + _random_dna(100, 12)  # some guides fail (homopolymer)

    def test_include_low_quality(self):
        from analyses import crispr as C
        everything = C.design_guides(self.SEQUENCE, config=C.DesignConfig(include_low_quality=True))
        failing = [c for c in everything.candidates if not c.is_high_quality]
        self.assertTrue(failing)
        filtered = C.design_guides(self.SEQUENCE, config=C.DesignConfig(include_low_quality=False))
        self.assertTrue(all(c.is_high_quality for c in filtered.candidates))
        self.assertEqual(len(filtered.candidates), len(everything.candidates) - len(failing))
        self.assertTrue(any("FAIL flags were excluded" in w for w in filtered.warnings))

    def test_rank_by_gc(self):
        from analyses import crispr as C
        with_gc = C.design_guides(self.SEQUENCE, config=C.DesignConfig(rank_by_gc=True))
        without_gc = C.design_guides(self.SEQUENCE, config=C.DesignConfig(rank_by_gc=False))
        scores = lambda r: {c.guide_start: c.rank_score for c in r.candidates}
        self.assertNotEqual(scores(with_gc), scores(without_gc))
        # GC is still reported as a component.
        self.assertIn("GC", without_gc.candidates[0].rank_components)

    def test_configured_editing_window_and_editor_are_used(self):
        from analyses import crispr as C
        mutation = C.MutationTarget("chr1", 150, self.SEQUENCE[149], "T" if self.SEQUENCE[149] != "T" else "C")
        default = C.config_editing_window(C.DesignConfig())
        self.assertEqual(default, tuple(C.BASE_EDITOR_DEFINITIONS["CBE"]["window"]))  # behavior unchanged
        from unittest import mock
        real = C.analyze_base_editing
        with mock.patch.object(C, "analyze_base_editing", wraps=real) as spy:
            C.find_base_editing_candidates(
                self.SEQUENCE, mutation, config=C.DesignConfig(editing_window_start=1, editing_window_end=2,
                                                               base_editor="ABE"))
        self.assertEqual((spy.call_args.kwargs["editor"], spy.call_args.kwargs["window"]), ("ABE", (1, 2)))
        analyses = real(C.design_guides(self.SEQUENCE).candidates, mutation, editor="ABE", window=(1, 2))
        self.assertTrue(analyses)
        self.assertTrue(all(item.editing_window == (1, 2) for item in analyses))
        with self.assertRaises(C.CRISPRAnalysisError):
            C.config_editing_window(C.DesignConfig(editing_window_start=5, editing_window_end=4))

    def test_relative_exports_go_to_the_project_or_the_configured_directory(self):
        import os
        import tempfile
        from pathlib import Path
        from unittest import mock
        from analyses import crispr as C
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, {"SYNAPSE_EXPORT_DIR": tmp}):
            self.assertEqual(C._default_export_directory(), Path(tmp) / "crispr")
        with mock.patch.dict(os.environ, {"SYNAPSE_EXPORT_DIR": ""}):
            self.assertEqual(C._default_export_directory(), PROJECT_ROOT / "exports" / "crispr")


# ----------------------------------------------------------------------
# D19: the pipeline's ORF runner takes the menu's search settings
# ----------------------------------------------------------------------
class PipelineOrfRunnerTests(unittest.TestCase):
    SEQUENCE = _random_dna(3000, 5)

    def test_run_orf_reproduces_api_results_for_every_setting(self):
        from core import api
        from core.pipeline import SynapsePipeline
        settings = {"start_codon_policy": "table_starts", "nested": True, "ranking": "heuristic",
                    "genetic_code": 11, "min_length": 60}
        pipeline = SynapsePipeline(sequence=self.SEQUENCE, header="x")
        result = pipeline.run_orf(**settings)
        expected = api.find_orfs(self.SEQUENCE, **settings).data
        self.assertEqual(result["parameters"], expected["parameters"])
        self.assertEqual([o["sequence"] for o in result["orfs"]], [o["sequence"] for o in expected["orfs"]])
        default = SynapsePipeline(sequence=self.SEQUENCE, header="x").run_orf()
        self.assertNotEqual(default["parameters"], result["parameters"])

    def test_completed_has_a_fixed_shape(self):
        from core.pipeline import STATUS_MODULES, SynapsePipeline
        pipeline = SynapsePipeline(sequence=self.SEQUENCE, header="x")
        self.assertEqual(tuple(pipeline.context.status()), STATUS_MODULES)
        pipeline.store("blast", {"hits": []})
        status = pipeline.context.status()
        self.assertEqual(tuple(status), STATUS_MODULES)
        self.assertTrue(status["blast"])
        pipeline.context.clear_results()
        self.assertEqual(tuple(pipeline.context.status()), STATUS_MODULES)


# ----------------------------------------------------------------------
# D15: the unimplemented planner is labeled as such
# ----------------------------------------------------------------------
class PlannerPlaceholderTests(unittest.TestCase):
    def test_menu_says_not_implemented(self):
        import os
        import subprocess
        import sys
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, SYNAPSE_LOG_DIR=tmp, PYTHONIOENCODING="utf-8")
            completed = subprocess.run([sys.executable, "main.py"], cwd=PROJECT_ROOT, input="\n1\n11\n\n0\n",
                                       capture_output=True, text=True, encoding="utf-8", timeout=120, env=env)
        output = completed.stdout + completed.stderr
        self.assertEqual(completed.returncode, 0, output[-2000:])
        self.assertIn("[11] AI Experiment Planner (not implemented)", output)
        self.assertIn("AI Experiment Planner: not implemented. Nothing is run or stored.", output)
        self.assertNotIn("coming soon", output)


# ----------------------------------------------------------------------
# S15, S16: variant module
# ----------------------------------------------------------------------
class VariantModuleFindingTests(unittest.TestCase):
    def test_ref_mismatch_flag_does_not_leak_into_the_input(self):
        from analyses import variant_analysis as V
        original = V.Variant("chr1", 1, "TTT", "T")
        original.annotation["validation_flags"] = ["PRE_EXISTING"]
        normalized = V.left_normalize_variant(original, "GGGAAAAAGGG")
        self.assertEqual(original.annotation["validation_flags"], ["PRE_EXISTING"])  # was + REF_MISMATCH
        self.assertEqual(normalized.annotation["validation_flags"], ["PRE_EXISTING", "REF_MISMATCH"])

    def test_mutation_bridge_runs(self):
        from analyses import variant_analysis as V
        variant = V.annotate_from_mutation_analysis(
            V.Variant("c", 4, "A", "G"), reference_sequence="ATGAAATAA", alternative_sequence="ATGGAATAA")
        self.assertEqual((variant.variant_type, variant.consequence, variant.coding_change),
                         ("SNV", "Missense", "AAA>GAA"))


# ----------------------------------------------------------------------
# S17: the cloning digest cuts where the enzyme cuts
# ----------------------------------------------------------------------
class CloningDigestCutTests(unittest.TestCase):
    def test_audit_example(self):
        from analyses import cloning as CL
        sequence = _random_dna(50, 21).replace("GAATTC", "GAATTA") + "GAATTC" + "A" * 50 + "CTGCAG" + "T" * 50
        digest = CL._selected_digest(sequence, ["EcoRI", "PstI"], "linear")
        # EcoRI G^AATTC at 51, PstI CTGCA^G at 107: cuts after 51 and 111 (was 51 and 107).
        self.assertEqual(sorted(f["length"] for f in digest["fragments"]), [51, 51, 60])

    def test_agrees_with_the_restriction_module_for_every_enzyme_pair(self):
        import itertools
        from analyses import cloning as CL, restriction as R
        from data.restriction_enzymes import RESTRICTION_ENZYMES
        sequence = _random_dna(1500, 9)
        for topology in ("linear", "circular"):
            for pair in itertools.combinations(sorted(RESTRICTION_ENZYMES), 2):
                with self.subTest(topology=topology, pair=pair):
                    cloning = CL._selected_digest(sequence, list(pair), topology)
                    restriction = R.perform_digest(sequence, list(pair), topology=topology)
                    self.assertEqual(sorted(f["length"] for f in cloning["fragments"]),
                                     sorted(f["length"] for f in restriction["fragments"]))


# ----------------------------------------------------------------------
# S19: start-codon classification and scores follow the genetic code
# ----------------------------------------------------------------------
class StartCodonClassificationTests(unittest.TestCase):
    def test_classification_uses_the_table(self):
        from analyses import orf as O
        self.assertEqual(O.classify_start_codon("ATG"), "Canonical")
        self.assertEqual(O.classify_start_codon("CTG"), "Alternative")      # NCBI table 1 start; was "None"
        self.assertEqual(O.classify_start_codon("ATA", 2), "Alternative")   # vertebrate mitochondrial
        self.assertEqual(O.classify_start_codon("ATA"), "None")
        self.assertEqual(O.classify_start_codon("AAA"), "None")

    def test_ctg_orf_scores_like_gtg_orf(self):
        from analyses import orf as O
        body = "GCT" * 40 + "TAA"
        ctg = next(o for o in O.find_orfs("CC" + "CTG" + body, min_length=60, start_codon_policy="table_starts")
                   if o["start_codon"] == "CTG")
        gtg = next(o for o in O.find_orfs("CC" + "GTG" + body, min_length=60,
                                          start_codon_policy="atg_and_alternative")
                   if o["start_codon"] == "GTG")
        self.assertEqual(O.calculate_cds_quality(ctg)["score"], O.calculate_cds_quality(gtg)["score"])
        self.assertIn("Alternative start codon", O.calculate_cds_quality(ctg)["reasons"])


# ----------------------------------------------------------------------
# S23: the stored guide seed is the PAM-proximal end of the guide
# ----------------------------------------------------------------------
class CrisprSeedOrientationTests(unittest.TestCase):
    SEQUENCE = _random_dna(400, 4)

    def test_seed_is_pam_proximal_for_every_system(self):
        from analyses import crispr as C
        for name in ("SpCas9", "SpCas9-NG", "SaCas9", "Cas12a"):
            system = C.get_crispr_system(name)
            result = C.design_guides(self.SEQUENCE, system=name)
            self.assertTrue(result.candidates)
            for candidate in result.candidates:
                with self.subTest(system=name, candidate=candidate.candidate_id):
                    seed = candidate.seed_sequence
                    self.assertEqual(len(seed), min(C.DEFAULT_SEED_LENGTH, len(candidate.guide)))
                    # The seed must lie directly against the PAM in target_sequence.
                    target, pam = candidate.target_sequence, candidate.pam
                    if system.pam_side == "3prime":
                        self.assertEqual(seed, candidate.guide[-len(seed):])
                        self.assertTrue(target.endswith(pam))
                        self.assertTrue(target[:-len(pam)].endswith(seed))
                    else:
                        self.assertEqual(seed, candidate.guide[:len(seed)])
                        self.assertTrue(target.startswith(pam))
                        self.assertTrue(target[len(pam):].startswith(seed))

    def test_cas12a_seed_is_the_five_prime_end_not_the_three_prime_end(self):
        """The defect stored the 3' end for a 5'-PAM system."""
        from analyses import crispr as C
        candidates = C.design_guides(self.SEQUENCE, system="Cas12a").candidates
        length = C.DEFAULT_SEED_LENGTH
        distinguishing = [c for c in candidates if c.guide[:length] != c.guide[-length:]]
        self.assertTrue(distinguishing, "fixture produced no guide with distinguishable ends")
        for candidate in distinguishing:
            with self.subTest(candidate=candidate.candidate_id):
                self.assertEqual(candidate.seed_sequence, candidate.guide[:length])
                self.assertNotEqual(candidate.seed_sequence, candidate.guide[-length:])

    def test_stored_seed_matches_the_specificity_screen(self):
        from analyses import crispr as C
        for name in ("SpCas9", "Cas12a"):
            with self.subTest(system=name):
                system = C.get_crispr_system(name)
                candidate = C.design_guides(self.SEQUENCE, system=name).candidates[0]
                self.assertEqual(
                    C.pam_proximal_seed(candidate.guide, system, C.DEFAULT_SEED_LENGTH),
                    candidate.seed_sequence)
                hits = C.search_local_off_targets_detailed(
                    candidate, {"self": self.SEQUENCE}, max_mismatches=0)
                self.assertTrue(hits)
                self.assertEqual([h["seed_mismatches"] for h in hits], [0] * len(hits))

    def test_seed_length_is_bounded(self):
        from analyses import crispr as C
        system = C.get_crispr_system("SpCas9")
        self.assertEqual(C.pam_proximal_seed("ACGTACGT", system, 0), "")
        self.assertEqual(C.pam_proximal_seed("ACGT", system, 99), "ACGT")


# ----------------------------------------------------------------------
# S3: edit_in_rtt reports actual coverage of the edit site
# ----------------------------------------------------------------------
class PrimeEditingRttTests(unittest.TestCase):
    SEQUENCE = _random_dna(400, 4)
    EDIT_POSITION = 200

    def _mutation(self):
        from analyses import crispr as C
        reference = self.SEQUENCE[self.EDIT_POSITION - 1]
        return C.MutationTarget("chr1", self.EDIT_POSITION, reference,
                                "T" if reference != "T" else "G")

    def _candidates(self):
        from analyses import crispr as C
        mutation = self._mutation()
        guides = C.design_guides(self.SEQUENCE, mutations=[mutation])
        return mutation, C.design_prime_editing_candidates(
            self.SEQUENCE, mutation, guides.candidates,
            pbs_lengths=(13,), rtt_lengths=(16,))

    def test_rtt_actually_covers_the_edit_site(self):
        from analyses import crispr as C
        mutation, items = self._candidates()
        self.assertTrue(items)
        for item in items:
            with self.subTest(candidate=item.candidate_id):
                self.assertTrue(item.edit_in_rtt)
                # The RTT is built from the edit site outwards, so its first
                # base is the edit-site base in the guide's orientation.
                expected = (mutation.reference if item.strand == "+"
                            else C.reverse_complement(mutation.reference))
                self.assertEqual(item.rtt_sequence[0], expected)

    def test_notes_state_that_the_rtt_does_not_encode_the_edit(self):
        _, items = self._candidates()
        self.assertTrue(items)
        for item in items:
            with self.subTest(candidate=item.candidate_id):
                self.assertTrue(any("does not encode the alternative allele" in note
                                    for note in item.sequence_level_notes))

    def test_an_rtt_that_misses_the_edit_site_is_rejected(self):
        """A length check alone cannot detect this; the coordinate check can."""
        from analyses import crispr as C
        mutation, items = self._candidates()
        candidate = C.design_guides(self.SEQUENCE, mutations=[mutation]).candidates[0]
        # Anchor the edit far outside the supplied sequence: a full-length RTT
        # can still be cut from the sequence, but it cannot cover the edit.
        distant = C.MutationTarget("chr1", len(self.SEQUENCE) + 500,
                                   mutation.reference, mutation.alternative)
        item = C._prime_edit_for_candidate(self.SEQUENCE, candidate, distant, 13, 16)
        self.assertFalse(item.edit_in_rtt)


# ----------------------------------------------------------------------
# S4, S5: documented thresholds and a precise poly-T warning
# ----------------------------------------------------------------------
class CrisprThresholdDocumentationTests(unittest.TestCase):
    def test_poly_t_warning_names_the_pol_iii_context(self):
        from analyses import crispr as C
        self.assertEqual(C.POL3_TERMINATOR_T_RUN, 4)
        metrics = C.characterize_guide("ACGTTTTACGTACGTACGT")
        self.assertTrue(metrics.poly_t)
        self.assertIn("WARN: poly-T run (possible U6/Pol III terminator)",
                      C.quality_flags(metrics, C.DesignConfig()))
        # Below the screening threshold.
        self.assertFalse(C.characterize_guide("ACGTTTACGTACGTACGTA").poly_t)

    def test_methods_documents_the_screening_thresholds(self):
        methods = (PROJECT_ROOT / "docs" / "METHODS.md").read_text(encoding="utf-8")
        for fragment in ("PMID 24476820", "PMID 27096365", "PMID 29160308", "PMID 33811918"):
            with self.subTest(fragment=fragment):
                self.assertIn(fragment, methods)


# ----------------------------------------------------------------------
# S21, S10: the GenBank reader reports what it could not read
# ----------------------------------------------------------------------
class GenBankReaderNoteTests(unittest.TestCase):
    RECORD = (
        "LOCUS       ODD                       20 bp    DNA     linear   SYN 31-OCT-2026\n"
        "DEFINITION  Synthetic record exercising unparsed forms.\n"
        "ACCESSION   ODD\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..20\n"
        "     misc_feature    order(1..5,10..15)\n"
        "     variation       AB123.1:3..8\n"
        "ORIGIN\n"
        "        1 acgtacgt-a cgtacgtacg\n"
        "//\n"
    )

    def _load(self, text):
        import tempfile
        from pathlib import Path
        from fileio.genbank import load_genbank
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "record.gb"
            path.write_text(text)
            return load_genbank(str(path))

    def test_symbols_dropped_from_origin_are_reported(self):
        [record] = self._load(self.RECORD)
        self.assertEqual(record["sequence"], "ACGTACGTACGTACGTACG")  # the '-' was dropped
        self.assertTrue(any("not sequence letters" in note and "'-'" in note
                            for note in record["notes"]), record["notes"])

    def test_unparsed_feature_locations_are_reported(self):
        [record] = self._load(self.RECORD)
        note = next(n for n in record["notes"] if "does not parse" in n)
        self.assertIn("2 feature location(s)", note)
        self.assertIn("misc_feature", note)
        self.assertIn("variation", note)
        unparsed = sorted(f["key"] for f in record["features"]
                          if not f["location_parsed"]["parsed"])
        self.assertEqual(unparsed, ["misc_feature", "variation"])

    def test_a_clean_record_produces_no_such_notes(self):
        from fileio.genbank import load_genbank
        [record] = load_genbank(str(PROJECT_ROOT / "tests" / "fixtures" / "synthetic_circular.gb"))
        self.assertEqual(record["notes"], [])
        self.assertTrue(all(f["location_parsed"]["parsed"] for f in record["features"]))


# ----------------------------------------------------------------------
# S8: a failing to_dict() is no longer swallowed silently
# ----------------------------------------------------------------------
class ReportSerializationTests(unittest.TestCase):
    def test_failing_to_dict_is_logged(self):
        from core.report import to_jsonable

        class Hostile:
            def to_dict(self):
                raise RuntimeError("boom")

        with self.assertLogs("synapse.report", level="WARNING") as captured:
            result = to_jsonable(Hostile())
        self.assertIsInstance(result, str)
        self.assertTrue(any("to_dict() failed for Hostile" in line for line in captured.output),
                        captured.output)

    def test_ordinary_payloads_are_unchanged(self):
        from core.report import to_jsonable
        self.assertEqual(to_jsonable({"a": (1, 2), "b": {3}}), {"a": [1, 2], "b": [3]})
        self.assertEqual(to_jsonable(float("nan")), "NaN")


# ----------------------------------------------------------------------
# S22: the supplied variant reference must be DNA
# ----------------------------------------------------------------------
class VariantReferenceAlphabetTests(unittest.TestCase):
    def test_a_reference_in_another_alphabet_is_refused(self):
        from analyses import variant_analysis as V
        variants = [V.Variant("chr1", 2, "C", "T")]
        for reference in ("MKTAYIAKQRQISFVK", "ACGU", "ACGT-ACGT"):
            with self.subTest(reference=reference):
                with self.assertRaisesRegex(V.VariantAnalysisError, "is not DNA"):
                    V.check_reference_alleles(variants, {"chr1": reference})
                with self.assertRaisesRegex(V.VariantAnalysisError, "is not DNA"):
                    V.normalize_variants(variants, contig_references={"chr1": reference})
        with self.assertRaisesRegex(V.VariantAnalysisError, "is not DNA"):
            V.normalize_variants(variants, reference_sequence="MKTAYIAKQRQISFVK")

    def test_ambiguity_codes_in_a_reference_are_accepted(self):
        from analyses import variant_analysis as V
        reference = "ACGTNRYACGT"  # real assemblies contain N and other IUPAC codes
        self.assertEqual(V.validate_reference_sequence(reference), reference)
        report = V.check_reference_alleles([V.Variant("chr1", 1, "A", "G")], {"chr1": reference})
        self.assertEqual(report["checked"], 1)


# ----------------------------------------------------------------------
# S1, S2: minus-strand ORF context and longest-per-stop collapsing
# ----------------------------------------------------------------------
class MinusStrandOrfTests(unittest.TestCase):
    """Both findings were speculative; these tests record the actual behavior."""

    def test_kozak_context_of_a_minus_strand_orf_is_the_minus_strand_context(self):
        from analyses import orf as O
        from core.seqtools import reverse_complement
        # Build the intended minus-strand transcript, then place its reverse
        # complement in the parent, so the parent carries it on the minus strand.
        minus = "TTTT" + "GCCACC" + "ATG" + "GCT" * 21 + "TAA" + "TTTT"
        parent = reverse_complement(minus)
        orf = next(o for o in O.find_orfs(parent, min_length=60) if o["strand"] == "-")
        kozak = orf["start_context"]["kozak"]
        self.assertTrue(kozak["available"])
        self.assertEqual(kozak["context"], "GCCACCATGG")  # read on the minus strand
        self.assertEqual(kozak["classification"], "Strong")
        # No coordinate is reported, so none can be mapped to the wrong frame.
        self.assertEqual(sorted(kozak), ["available", "classification", "context", "score"])

    def test_longest_orf_per_stop_collapses_on_the_minus_strand(self):
        from analyses import orf as O
        from core.seqtools import reverse_complement
        minus = "ATG" + "GCT" * 10 + "ATG" + "GCT" * 10 + "TAA"
        parent = "CCCC" + reverse_complement(minus) + "CCCC"
        nested = [o for o in O.find_orfs(parent, min_length=30, nested=True)
                  if o["strand"] == "-"]
        collapsed = [o for o in O.find_orfs(parent, min_length=30, nested=False)
                     if o["strand"] == "-"]
        self.assertEqual(sorted(o["length"] for o in nested), [36, 69])
        # Both share one reading frame and one 3' end, so only the longest stays.
        self.assertEqual({o["frame"] for o in nested}, {-2})
        self.assertEqual({o["start"] for o in nested}, {4})
        self.assertEqual([o["length"] for o in collapsed], [69])


# ----------------------------------------------------------------------
# S9: the alignment note states the actual cost
# ----------------------------------------------------------------------
class AlignmentCostNoteTests(unittest.TestCase):
    def test_note_quantifies_the_memory_cost(self):
        from analyses.mutation import MAX_RECOMMENDED_ALIGNMENT_LENGTH, compare_sequences
        length = MAX_RECOMMENDED_ALIGNMENT_LENGTH + 10
        reference = "ACGT" * (length // 4)
        variant = reference[:100] + "T" + reference[101:]
        note = next(n for n in compare_sequences(reference, variant).notes
                    if "aligned globally" in n)
        self.assertIn("bytes per cell", note)
        self.assertIn("4.5 GB", note)

    def test_short_sequences_get_no_such_note(self):
        from analyses.mutation import compare_sequences
        reference = "ACGT" * 25
        result = compare_sequences(reference, reference[:50] + "T" + reference[51:])
        self.assertFalse(any("aligned globally" in n for n in result.notes))


# ----------------------------------------------------------------------
# P1-9: the RNA codon view uses the selected ORF's genetic code
# ----------------------------------------------------------------------
class RnaGeneticCodeTests(unittest.TestCase):
    # Table 2 reads UGA as Trp, table 1 as a stop.
    RNA = "AUG" + "UGA" * 5 + "AAA" + "UAA"
    DNA = RNA.replace("U", "T")

    def test_frame_translation_follows_the_table(self):
        from analyses import rna as R
        standard = R.translate_rna_frame(self.RNA, 0)
        mitochondrial = R.translate_rna_frame(self.RNA, 0, genetic_code=2)
        self.assertEqual(standard["protein"], "M")              # stops at the first UGA
        self.assertEqual(mitochondrial["protein"], "MWWWWWK")   # UGA is Trp
        self.assertEqual((standard["genetic_code"], mitochondrial["genetic_code"]), (1, 2))

    def test_three_frame_view_records_the_table(self):
        from analyses import rna as R
        view = R.analyze_rna_three_frames(self.RNA, genetic_code=2)
        self.assertEqual(view["genetic_code"], 2)
        self.assertTrue(all(f["genetic_code"] == 2 for f in view["frames"]))
        self.assertEqual(view["frames"][0]["protein"], "MWWWWWK")

    def test_a_selected_orf_supplies_its_table(self):
        from analyses import rna as R
        orf = {"sequence": self.DNA, "genetic_code": 2}
        self.assertEqual(R.orf_genetic_code(orf), 2)
        self.assertEqual(R.orf_genetic_code({"sequence": self.DNA}), 1)
        self.assertEqual(R.orf_genetic_code(None), 1)
        analysis = R.analyze_rna_from_pipeline(self.DNA, selected_orf=orf)
        self.assertEqual(analysis["genetic_code"], 2)
        self.assertEqual(analysis["translation_frames"]["frames"][0]["protein"], "MWWWWWK")

    def test_without_an_orf_the_default_table_is_used(self):
        from analyses import rna as R
        analysis = R.analyze_rna_from_pipeline(self.DNA)
        self.assertEqual(analysis["genetic_code"], 1)
        self.assertEqual(analysis["translation_frames"]["frames"][0]["protein"], "M")

    def test_api_and_orf_translation_agree(self):
        from core import api
        from analyses.orf import translate_orf_sequence
        orf = {"sequence": self.DNA, "genetic_code": 2}
        view = api.analyze_rna(self.DNA, genetic_code=2).data["translation_frames"]
        # The +1 frame preview must match what the ORF itself translates to.
        self.assertEqual(view["frames"][0]["protein"], translate_orf_sequence(orf))


# ----------------------------------------------------------------------
# S18: the cloning planner no longer repeats identical site scans
# ----------------------------------------------------------------------
class CloningPlannerScanReuseTests(unittest.TestCase):
    VECTOR = _random_dna(1300, 5)
    INSERT = _random_dna(300, 6)

    def test_each_distinct_scan_runs_once(self):
        from analyses import cloning as C, restriction as R
        from data.restriction_enzymes import RESTRICTION_ENZYMES
        enzymes = sorted(RESTRICTION_ENZYMES)[:8]
        R._scan_recognition_sites.cache_clear()
        C.automatic_cloning_planner(self.VECTOR, self.INSERT, enzyme_names=enzymes)
        info = R._scan_recognition_sites.cache_info()
        # Two sequences x 8 enzymes x one topology each.
        self.assertEqual(info.misses, 2 * len(enzymes))
        # The planner asked far more often than that; the rest was reused.
        self.assertGreater(info.hits, 100 * info.misses)

    def test_a_cached_result_cannot_be_corrupted_by_a_caller(self):
        from analyses import restriction as R
        sequence = "GAATTC" + "A" * 20 + "GAATTC"
        first = R.find_recognition_site_matches(sequence, "GAATTC")
        first[0]["status"] = "TAMPERED"
        first.append({"bogus": True})
        second = R.find_recognition_site_matches(sequence, "GAATTC")
        self.assertEqual([m["status"] for m in second], ["CONFIRMED", "CONFIRMED"])
        self.assertEqual(len(second), 2)

    def test_caching_does_not_change_results(self):
        from analyses import restriction as R
        sequence = "GAATTCNNGAATTC" + _random_dna(200, 7)
        for topology in ("linear", "circular"):
            with self.subTest(topology=topology):
                R._scan_recognition_sites.cache_clear()
                cold = R.find_recognition_site_matches(sequence, "GAATTC", topology=topology)
                warm = R.find_recognition_site_matches(sequence, "GAATTC", topology=topology)
                self.assertEqual(cold, warm)
                self.assertGreater(R._scan_recognition_sites.cache_info().hits, 0)


# ----------------------------------------------------------------------
# S14: the API covers cloning, variant filtering, mutation reports,
#      BLAST and literature, without weakening the privacy gate
# ----------------------------------------------------------------------
class ApiCompletenessTests(unittest.TestCase):
    # Repeating fillers that cannot spell GAATTC or CTGCAG, so each site below
    # is the only one of its kind and the planner has a unique cut to work with.
    VECTOR = "ACAC" * 75 + "GAATTC" + "ATAT" * 10 + "CTGCAG" + "AGAG" * 75
    INSERT = "GAATTC" + "ACGA" * 62 + "CTGCAG"

    def test_every_analysis_module_has_an_entry_point(self):
        from core import api
        for name in ("plan_cloning", "filter_variants", "mutation_report",
                     "blast_search", "literature_search"):
            with self.subTest(entry=name):
                self.assertTrue(callable(getattr(api, name)))
                self.assertIn(name, api.__all__)

    def test_plan_cloning_finds_a_valid_strategy(self):
        from core import api
        result = api.plan_cloning(self.VECTOR, self.INSERT,
                                  enzymes=["EcoRI", "PstI", "BamHI", "HindIII"])
        self.assertEqual(result.module, "cloning")
        self.assertEqual(result.data["valid_count"], 1)
        self.assertEqual(result.data["best"]["mode"], "two-enzyme")
        self.assertEqual(result.parameters["vector_topology"], "linear")

    def test_plan_cloning_refuses_a_non_dna_sequence(self):
        from core import api
        from core.alphabets import SequenceTypeError
        with self.assertRaisesRegex(SequenceTypeError, "needs a DNA sequence"):
            api.plan_cloning(self.VECTOR, self.INSERT.replace("T", "U"))

    def test_filter_variants(self):
        from core import api
        from analyses import variant_analysis as V
        variants = [V.Variant("chr1", 10, "A", "G", quality=50.0),
                    V.Variant("chr1", 20, "C", "T", quality=5.0)]
        result = api.filter_variants(variants, minimum_quality=20.0)
        self.assertEqual(result.module, "variant_filter")
        self.assertEqual((len(result.data.passed), len(result.data.failed)), (1, 1))
        self.assertEqual(result.parameters["input_variants"], 2)
        with self.assertRaisesRegex(ValueError, "not both"):
            api.filter_variants(variants, V.VariantFilter(), minimum_quality=1.0)

    def test_mutation_report_renders_text(self):
        from core import api
        result = api.mutation_report("ATGAAATTTGGGTAA", "ATGAAATTCGGGTAA")
        self.assertEqual(result.module, "mutation_report")
        self.assertIn("MUTATION ANALYSIS REPORT", result.data["report_text"])
        self.assertTrue(result.limitations)

    def test_remote_blast_requires_explicit_confirmation(self):
        from core import api
        with self.assertRaisesRegex(ValueError, "send_to_ncbi=True"):
            api.blast_search("ACGTACGTACGTACGTACGT", engine="remote")

    def test_literature_search_requires_explicit_confirmation(self):
        from core import api
        with self.assertRaisesRegex(ValueError, "send_to_providers=True"):
            api.literature_search("CRISPR off-target")

    @unittest.skipUnless(_blast_plus_available("blastn", "blastdbcmd"),
                         "BLAST+ (blastn, blastdbcmd) not installed")
    def test_local_blast_search_sends_nothing(self):
        """External tool execution: the local engine needs no confirmation."""
        from core import api
        lines = (PROJECT_ROOT / "data" / "local_test" / "query.fasta").read_text().splitlines()
        query = "".join(line.strip() for line in lines if not line.startswith(">"))
        result = api.blast_search(query, engine="local", program="blastn",
                                  local_database_path=str(LOCAL_TEST_DB))
        self.assertEqual(result.module, "blast")
        self.assertEqual([hit.hit_id for hit in result.data.hits],
                         ["Test_Gene_001", "Test_Gene_002"])
        self.assertEqual(result.warnings, [])


# ----------------------------------------------------------------------
# S6 / P1-7: resuming an NCBI search from its request ID
# ----------------------------------------------------------------------
class BlastResumeTests(unittest.TestCase):
    """Driven by a stubbed client: no request leaves this machine."""

    RID = "ABCD1234567"

    def _client(self, status_sequence=("READY",), xml=None):
        from analyses import blast as B
        client = B.NCBIBlastClient()
        statuses = list(status_sequence)
        calls = {"status": 0, "fetch": 0, "submit": 0}

        def fake_status(rid):
            calls["status"] += 1
            return (statuses.pop(0) if statuses else "READY"), "stub"

        def fake_fetch(rid):
            calls["fetch"] += 1
            return xml if xml is not None else B.SAMPLE_XML

        def fake_submit(*args, **kwargs):
            calls["submit"] += 1
            raise AssertionError("resume() must not submit a new search")

        client.status = fake_status          # type: ignore[method-assign]
        client.fetch_xml = fake_fetch        # type: ignore[method-assign]
        client.submit = fake_submit          # type: ignore[method-assign]
        client.min_poll_interval = 0
        return client, calls

    def test_resume_fetches_without_resubmitting(self):
        client, calls = self._client()
        run = client.resume(self.RID, poll_interval=0)
        self.assertEqual(run.metadata.rid, self.RID)
        self.assertEqual(calls["submit"], 0)
        self.assertEqual(calls["fetch"], 1)
        self.assertTrue(any("resumed from its request ID" in w
                            for w in run.metadata.warnings))

    # NCBI output carries the query definition, length, program and database.
    XML_WITH_QUERY = (
        '<?xml version="1.0"?>\n<BlastOutput>\n'
        "  <BlastOutput_program>blastp</BlastOutput_program>\n"
        "  <BlastOutput_db>swissprot</BlastOutput_db>\n"
        "  <BlastOutput_query-def>SYN_TEST some description</BlastOutput_query-def>\n"
        "  <BlastOutput_query-len>137</BlastOutput_query-len>\n"
        "  <BlastOutput_iterations></BlastOutput_iterations>\n</BlastOutput>\n"
    )

    def test_query_and_parameters_are_recovered_from_the_xml(self):
        from analyses import blast as B
        query, parameters = B.describe_run_from_xml(self.XML_WITH_QUERY)
        self.assertEqual((parameters.program, parameters.database), ("blastp", "swissprot"))
        self.assertEqual(query.name, "SYN_TEST")
        self.assertEqual(query.sequence_type, "PROTEIN")
        self.assertEqual(query.sequence, "")          # never invented
        self.assertEqual((query.length, query.reported_length), (137, 137))
        self.assertTrue(any("not part of it" in note for note in query.notes))

    def test_a_missing_query_length_is_stated_not_guessed(self):
        from analyses import blast as B
        query, _ = B.describe_run_from_xml(B.SAMPLE_XML)  # carries no query-len
        self.assertEqual(query.length, 0)
        self.assertTrue(any("no query length" in note for note in query.notes))

    def test_a_supplied_query_is_preferred_over_the_xml(self):
        from analyses import blast as B
        client, _ = self._client()
        original = B.build_query("ACGTACGTACGTACGTACGT", name="Original")
        run = client.resume(self.RID, query=original,
                            parameters=B.BlastParameters(program="blastn", database="nt"),
                            submitted_at="2026-10-31T00:00:00+00:00", poll_interval=0)
        self.assertEqual(run.query.name, "Original")
        self.assertEqual(run.query.sequence, original.sequence)
        self.assertFalse(any("resumed from its request ID" in w
                             for w in run.metadata.warnings))

    def test_waiting_then_ready(self):
        client, calls = self._client(status_sequence=("WAITING", "WAITING", "READY"))
        client.resume(self.RID, poll_interval=0)
        self.assertEqual(calls["status"], 3)

    def test_an_expired_request_is_explained(self):
        from analyses import blast as B
        client, _ = self._client(status_sequence=("UNKNOWN",))
        with self.assertRaisesRegex(B.BlastAnalysisError, "may simply have expired"):
            client.resume(self.RID, poll_interval=0)

    def test_an_empty_request_id_is_refused(self):
        from analyses import blast as B
        client, _ = self._client()
        with self.assertRaisesRegex(B.BlastAnalysisError, "request ID is required"):
            client.resume("   ")

    def test_api_requires_confirmation_before_contacting_ncbi(self):
        from core import api
        with self.assertRaisesRegex(ValueError, "contact_ncbi=True"):
            api.blast_resume(self.RID)

    def test_the_menu_offers_the_entry(self):
        import inspect
        from analyses import blast as B
        source = inspect.getsource(B.blast_menu)
        self.assertIn("Resume Search by Request ID", source)
        self.assertIn('choice == "2"', source)


# ----------------------------------------------------------------------
# S13: importing main must not start the application
# ----------------------------------------------------------------------
class MainEntryPointTests(unittest.TestCase):
    def test_importing_main_does_not_run_the_session(self):
        import os
        import subprocess
        import sys
        import tempfile
        script = (
            "import main\n"
            "print('IMPORTED', callable(main.main))\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, SYNAPSE_LOG_DIR=tmp, PYTHONIOENCODING="utf-8")
            # No stdin: if the module still asked for input, this would fail.
            completed = subprocess.run(
                [sys.executable, "-c", script], cwd=PROJECT_ROOT, input="",
                capture_output=True, text=True, encoding="utf-8", timeout=120, env=env,
            )
        output = completed.stdout + completed.stderr
        self.assertEqual(completed.returncode, 0, output[-2000:])
        self.assertIn("IMPORTED True", completed.stdout)
        self.assertNotIn("SYNAPSE - LOAD SEQUENCE FILE", output)
        self.assertNotIn("Traceback", output)

    def test_module_level_flow_is_behind_the_guard(self):
        import ast
        tree = ast.parse((PROJECT_ROOT / "main.py").read_text(encoding="utf-8"))
        loops = [n for n in tree.body if isinstance(n, (ast.While, ast.For))]
        self.assertEqual(loops, [], "main.py must not loop at module level")
        guards = [
            n for n in tree.body
            if isinstance(n, ast.If) and ast.unparse(n.test) == "__name__ == '__main__'"
        ]
        self.assertEqual(len(guards), 1)
        self.assertIn("main()", ast.unparse(guards[0]))


# ----------------------------------------------------------------------
# S24: what the BLAST XML parser is and is not exposed to
# ----------------------------------------------------------------------
class XmlParsingExposureTests(unittest.TestCase):
    """Turns the S24 assumption into measurements that hold across versions."""

    HEAD = '<?xml version="1.0"?>\n'
    BODY = "<BlastOutput><BlastOutput_program>{}</BlastOutput_program></BlastOutput>"

    def test_external_entities_are_not_resolved(self):
        """XXE: the stdlib parser never fetches an external entity."""
        import xml.etree.ElementTree as ET
        from analyses.blast import BlastAnalysisError, describe_run_from_xml
        payload = (self.HEAD
                   + '<!DOCTYPE r [<!ENTITY xxe SYSTEM "file:///C:/Windows/win.ini">]>\n'
                   + self.BODY.format("&xxe;"))
        with self.assertRaises(ET.ParseError) as caught:
            ET.fromstring(payload)
        self.assertIn("undefined entity", str(caught.exception))
        # The module turns that into its own error rather than leaking a parse error.
        with self.assertRaises(BlastAnalysisError):
            describe_run_from_xml(payload)

    def test_an_external_dtd_is_not_fetched(self):
        import xml.etree.ElementTree as ET
        payload = (self.HEAD
                   + '<!DOCTYPE r SYSTEM "http://example.invalid/evil.dtd">\n'
                   + "<BlastOutput/>")
        # Parses, and the unreachable host is never contacted.
        self.assertEqual(ET.fromstring(payload).tag, "BlastOutput")

    @staticmethod
    def _nested_entities(levels, fan=10):
        entities = ['<!ENTITY e0 "aaaaaaaaaa">']
        for index in range(1, levels):
            entities.append(f'<!ENTITY e{index} "{("&e%d;" % (index - 1)) * fan}">')
        return ("<?xml version=\"1.0\"?>\n<!DOCTYPE r [\n" + "\n".join(entities)
                + f"\n]>\n<BlastOutput><BlastOutput_program>&e{levels - 1};"
                  "</BlastOutput_program></BlastOutput>")

    def test_entity_expansion_is_bounded_by_expat(self):
        """Billion laughs is reachable but capped, when Expat is new enough."""
        import xml.etree.ElementTree as ET
        import xml.parsers.expat as expat
        if expat.version_info < (2, 4, 0):
            self.skipTest(f"Expat {expat.version_info} predates the amplification guard")
        with self.assertRaises(ET.ParseError) as caught:
            ET.fromstring(self._nested_entities(8))
        self.assertIn("amplification", str(caught.exception))

    def test_a_realistic_document_still_parses(self):
        from analyses import blast as B
        query, parameters = B.describe_run_from_xml(B.SAMPLE_XML)
        self.assertEqual(parameters.program, "blastn")
        self.assertIsNotNone(query)


# ----------------------------------------------------------------------
# D18: heuristic scores are not called "confidence" in data or exports
# ----------------------------------------------------------------------
def _confidence_keys(value, path="$"):
    """Every dictionary key below ``value`` that contains "confidence"."""
    found = []
    if isinstance(value, dict):
        for key, item in value.items():
            if "confidence" in str(key).lower():
                found.append(f"{path}.{key}")
            found.extend(_confidence_keys(item, f"{path}.{key}"))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_confidence_keys(item, f"{path}[{index}]"))
    return found


def _literal_keys_containing(word, directories):
    """String literals used as dict keys or subscripts that contain ``word``."""
    import ast
    hits = []
    for directory in directories:
        for path in sorted((PROJECT_ROOT / directory).glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                keys = []
                if isinstance(node, ast.Dict):
                    keys = [k for k in node.keys if k is not None]
                elif isinstance(node, ast.Subscript):
                    keys = [node.slice]
                for key in keys:
                    if (isinstance(key, ast.Constant) and isinstance(key.value, str)
                            and word in key.value.lower()):
                        hits.append(f"{directory}/{path.name}:{key.lineno} {key.value!r}")
    return hits


class HeuristicScoreVocabularyTests(unittest.TestCase):
    def test_orf_results_carry_no_confidence_key(self):
        from analyses import orf as O
        from tests.golden_scenarios import RANDOM_600
        result = O.analyze_orfs_enhanced(RANDOM_600)
        self.assertGreater(result["candidate_count"], 0)
        self.assertEqual(_confidence_keys(result), [])
        orf = result["orfs"][0]
        self.assertEqual(set(orf["heuristic_score"]), {"score", "classification"})
        self.assertEqual(_confidence_keys(O.analyze_selected_orf(orf)), [])

    def test_comparison_rows_and_summary_use_heuristic_names(self):
        from analyses import orf as O
        from tests.golden_scenarios import RANDOM_600
        orfs = O.analyze_orfs_enhanced(RANDOM_600)["orfs"]
        row = O.compare_orf_candidates(orfs, top_n=1)[0]
        self.assertEqual((row["heuristic_score"], row["heuristic_classification"]),
                         (orfs[0]["heuristic_score"]["score"],
                          orfs[0]["heuristic_score"]["classification"]))
        summary = O.calculate_orf_feature_summary(orfs)
        self.assertEqual(summary["high_heuristic_score"] + summary["moderate_heuristic_score"]
                         + summary["low_heuristic_score"], len(orfs))

    def test_former_function_name_still_works(self):
        from analyses import orf as O
        self.assertIs(O.calculate_orf_confidence, O.calculate_orf_heuristic_score)

    def test_signal_peptide_screen_reports_evidence(self):
        from analyses import protein as P
        signal = P.predict_signal_peptide("MKLLLLAAALLLLACSSA" + "K" * 20)
        self.assertIn(signal["evidence"], ("low", "moderate"))
        self.assertNotIn("confidence", signal)

    def test_no_confidence_literal_keys_in_source(self):
        directories = ("analyses", "core", "fileio", "ui", ".")
        self.assertEqual(_literal_keys_containing("confidence", directories), [])

    def test_sequence_classification_reports_an_evidence_level(self):
        from core.alphabets import EVIDENCE_LEVELS, classify_sequence
        from core.models import SequenceRecord
        result = classify_sequence("MKTAYIAKQR")
        self.assertIn(result.evidence, EVIDENCE_LEVELS)
        self.assertEqual(_confidence_keys(result.as_dict()), [])
        record = SequenceRecord("p", "MKTAYIAKQR")
        self.assertEqual(record.type_evidence, result.evidence)


# ----------------------------------------------------------------------
# Stage G: units and coordinate notation in printed output
# ----------------------------------------------------------------------
def _fstring_literal_parts(directories):
    """(location, text) for every literal part of every f-string."""
    import ast
    parts = []
    for directory in directories:
        for path in sorted((PROJECT_ROOT / directory).glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, ast.JoinedStr):
                    for value in node.values:
                        if isinstance(value, ast.Constant) and isinstance(value.value, str):
                            parts.append((f"{directory}/{path.name}:{node.lineno}", value.value))
    return parts


class OutputNotationTests(unittest.TestCase):
    DIRECTORIES = ("analyses", "core", "fileio", "ui", ".")

    def test_printed_units_use_one_notation(self):
        import re
        rules = {
            "space before percent sign": re.compile(r"^ %"),
            "no space before degree Celsius": re.compile(r"[^ ]°C"),
            "entropy unit other than cal/(mol·K)": re.compile(r"cal/mol/K|mol\*K"),
            "molar mass unit instead of Da": re.compile(r"\bg/mol\b"),
        }
        violations = [f"{where}: {name}: {text!r}"
                      for where, text in _fstring_literal_parts(self.DIRECTORIES)
                      for name, rule in rules.items() if rule.search(text)]
        self.assertEqual(violations, [])

    def test_blast_text_report_shows_e_values_with_three_significant_figures(self):
        import re
        from analyses import blast as B
        query = B.build_query("ATGCGTACGTACGTACGTAC" + "A" * 80, name="notation")
        run = B.parse_blast_xml(B.SAMPLE_XML, query=query, rid="NOTATION",
                                parameters=B.BlastParameters(),
                                submitted_at="2026-01-01T00:00:00+00:00")
        report = B.build_text_report(run)
        shown = re.findall(r"E-value:\s+(\S+)", report)
        self.assertTrue(shown)
        expected = [B._evalue_text(hit.best_evalue) for hit in B.rank_hits(run.hits)[:20]]
        self.assertEqual(shown, expected)
        self.assertEqual((B._evalue_text(1.234567e-50), B._evalue_text(None)), ("1.23e-50", "n/a"))

    def _pair(self):
        from analyses import primer as P
        from tests._support import deterministic_dna
        sequence = deterministic_dna(400)
        ranked = P.design_primers(sequence).ranked_pairs[0]
        return P, sequence, getattr(ranked, "pair", ranked)

    def test_primer_positions_are_shown_one_based(self):
        P, sequence, pair = self._pair()
        text = P.summarize_primer_pair(pair)
        forward = pair.forward
        self.assertEqual(sequence[forward.start:forward.end], forward.sequence)
        self.assertIn(f"Position:              {forward.start + 1}-{forward.end} (1-based)", text)

    def test_primer_target_input_is_one_based_inclusive(self):
        import builtins
        import contextlib
        import io
        from unittest import mock
        from analyses import primer as P
        sequence = "A" * 200
        cases = {("11", "60"): (10, 60), ("", ""): (0, 200), ("1", "200"): (0, 200),
                 ("0", "60"): None, ("61", "60"): None, ("1", "201"): None}
        for answers, expected in cases.items():
            with mock.patch.object(builtins, "input", side_effect=list(answers)), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(P._primer_ui_read_target(sequence), expected, answers)


# ----------------------------------------------------------------------
# Stage G (6.4): warnings a result carries reach the exported summary
# ----------------------------------------------------------------------
class ExportedWarningsTests(unittest.TestCase):
    def test_payload_warnings_reach_the_markdown_report(self):
        import logging
        import types
        from core import api
        from core.report import build_analysis_report, report_to_markdown
        from tests._support import deterministic_dna
        logging.disable(logging.CRITICAL)
        try:
            result = api.design_primers(deterministic_dna(600))
        finally:
            logging.disable(logging.NOTSET)
        payload_warnings = list(result.data.warnings)
        self.assertTrue(payload_warnings)
        self.assertTrue(set(payload_warnings) <= set(result.warnings))
        context = types.SimpleNamespace(results={"primer": [result]}, sequence="", header="")
        markdown = report_to_markdown(build_analysis_report(context))
        for warning in payload_warnings:
            self.assertIn(f"  - {warning}", markdown)

    def test_warnings_passed_and_carried_are_not_duplicated(self):
        from core.models import AnalysisResult
        payload = {"warnings": ["W1", "W2"]}
        result = AnalysisResult(module="x", data=payload, warnings=["W0", "W1"])
        self.assertEqual(result.warnings, ["W0", "W1", "W2"])
        self.assertEqual(AnalysisResult(module="x", data={"value": 1}).warnings, [])

    def test_restriction_results_record_parameters_and_enzyme_data_version(self):
        from analyses.restriction import analyze_restriction
        from core import api
        from core.pipeline import SynapsePipeline
        from data.restriction_enzymes import ENZYME_DATABASE_VERSION
        from tests._support import deterministic_dna
        sequence = deterministic_dna(600)
        pipeline = SynapsePipeline(sequence=sequence, header="x")
        pipeline.context.store_restriction_analysis(
            analyze_restriction(sequence, ["EcoRI", "BamHI"], "linear"))
        latest = pipeline.context.latest("restriction")
        assert latest is not None
        stored = latest.parameters
        self.assertEqual(stored, {"enzymes": ["EcoRI", "BamHI"], "topology": "linear",
                                  "enzyme_database_version": ENZYME_DATABASE_VERSION})
        digest = api.restriction_digest(sequence, ["EcoRI"])
        self.assertEqual(digest.parameters["enzyme_database_version"], ENZYME_DATABASE_VERSION)


# ----------------------------------------------------------------------
# Stage H (8.6, 9.6): one spelling, American English as in the identifiers
# ----------------------------------------------------------------------
def _python_sources():
    paths = [PROJECT_ROOT / "main.py"]
    for directory in ("analyses", "core", "data", "fileio", "ui", "tests"):
        paths.extend(sorted((PROJECT_ROOT / directory).glob("*.py")))
    return paths


class SpellingStandardTests(unittest.TestCase):
    BRITISH = ("behaviour", "analyse", "analysed", "labelled", "labelling", "modelled",
               "modelling", "normalise", "normalised", "normalisation", "organise",
               "organised", "organises", "serialise", "serialised", "serialisation",
               "recognise", "recognised", "unrecognised", "initialised", "optimisation",
               "prioritised", "prioritisation", "parameterised", "characterisation",
               "standardised", "standardisation", "summarise", "summarised", "colour",
               "centre")

    def test_text_and_documentation_use_american_spelling(self):
        import re
        pattern = re.compile(r"\b(" + "|".join(self.BRITISH) + r")\b", re.IGNORECASE)
        paths = [PROJECT_ROOT / "README.md", PROJECT_ROOT / "main.py"]
        for directory in ("docs", "analyses", "core", "fileio", "ui", "tests"):
            paths.extend(sorted((PROJECT_ROOT / directory).glob("*.md")))
            paths.extend(sorted((PROJECT_ROOT / directory).glob("*.py")))
        # docs/CHANGES.md must quote the original wording of an editorial
        # change (Master Prompt 8.7), so it is the one exempt document.
        exempt = {pathlib.Path(__file__).resolve(), (PROJECT_ROOT / "docs" / "CHANGES.md").resolve()}
        hits = []
        for path in paths:
            if path.resolve() in exempt:
                continue
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                hits.extend(f"{path.relative_to(PROJECT_ROOT).as_posix()}:{number} {m.group(0)}"
                            for m in pattern.finditer(line))
        self.assertEqual(hits, [])


# ----------------------------------------------------------------------
# Stage I (11.8): Python 3.10 compatibility without a 3.10 interpreter
# ----------------------------------------------------------------------
class Python310CompatibilityTests(unittest.TestCase):
    def test_source_parses_with_the_python_3_10_grammar(self):
        import ast
        failures = []
        for path in _python_sources():
            try:
                ast.parse(path.read_text(encoding="utf-8"), feature_version=(3, 10))
            except SyntaxError as exc:
                failures.append(f"{path.relative_to(PROJECT_ROOT).as_posix()}:{exc.lineno} {exc.msg}")
        self.assertEqual(failures, [])

    @unittest.skipUnless(hasattr(__import__("tokenize"), "FSTRING_START"),
                         "needs the Python 3.12+ tokenizer to see inside f-strings")
    def test_no_f_string_needs_python_3_12(self):
        """Before PEP 701 (3.12) an f-string field may not contain '#', a
        backslash, or the f-string's own quote. ``ast.parse`` with
        ``feature_version`` does not check this, and neither does pyright,
        which missed ``f"{'#':>3}"`` in the CRISPR module."""
        import io
        import tokenize
        # Looked up by name: the constants exist only from Python 3.12 on.
        fstring_start = getattr(tokenize, "FSTRING_START")
        fstring_end = getattr(tokenize, "FSTRING_END")
        hits = []
        for path in _python_sources():
            quotes = []
            source = path.read_text(encoding="utf-8")
            for token in tokenize.generate_tokens(io.StringIO(source).readline):
                if token.type == fstring_start:
                    if quotes and token.string.lstrip("rRfFbBuU")[0] == quotes[-1]:
                        hits.append((path, token.start[0], "nested f-string reuses the quote"))
                    quotes.append(token.string.lstrip("rRfFbBuU")[0])
                elif token.type == fstring_end:
                    quotes.pop()
                elif quotes and token.type == tokenize.COMMENT:
                    hits.append((path, token.start[0], "comment in an f-string field"))
                elif quotes and token.type == tokenize.STRING:
                    body = token.string.lstrip("rRbBuU")
                    if "#" in body or "\\" in body or body[0] == quotes[-1]:
                        hits.append((path, token.start[0], f"string {token.string!r} in an f-string field"))
        self.assertEqual([f"{p.relative_to(PROJECT_ROOT).as_posix()}:{line} {why}"
                          for p, line, why in hits], [])


class StyleGuideConsistencyTests(unittest.TestCase):
    def test_style_guide_names_only_tests_that_exist(self):
        import re
        guide = (PROJECT_ROOT / "docs" / "STYLE_GUIDE.md").read_text(encoding="utf-8")
        enforcement = guide.split("## 8. Enforcement", 1)[1]
        named = set(re.findall(r"`((?:[A-Z]\w+Tests)(?:\.test_\w+)?)`", enforcement))
        self.assertTrue(named)
        missing = []
        for name in sorted(named):
            class_name, _, method = name.partition(".")
            cls = globals().get(class_name)
            if cls is None or (method and not hasattr(cls, method)):
                missing.append(name)
        self.assertEqual(missing, [])
        self.assertTrue((PROJECT_ROOT / "tests" / "test_terminal_ui.py").exists())


if __name__ == "__main__":
    unittest.main()
