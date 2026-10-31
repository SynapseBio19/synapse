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


if __name__ == "__main__":
    unittest.main()
