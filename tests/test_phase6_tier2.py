"""Phase 6 Tier 2 tests: functional defects (docs/PHASE6_INITIAL_AUDIT.md)."""

from __future__ import annotations

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
# D8: modelled cut site, from cited conventions
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
        self.assertEqual(default, tuple(C.BASE_EDITOR_DEFINITIONS["CBE"]["window"]))  # behaviour unchanged
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
# D15: the unimplemented planner is labelled as such
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


if __name__ == "__main__":
    unittest.main()
