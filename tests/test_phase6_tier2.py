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


if __name__ == "__main__":
    unittest.main()
