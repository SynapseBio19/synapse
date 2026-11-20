"""The interval sweep in analyze_orf_relationships equals the all-pairs loop.

Until 2026-11-20 every ORF was compared with every other one, which took
about 8 minutes for the 32,406 ORFs of the E. coli genome
(docs/BENCHMARKS.md). The sweep only visits ORFs whose start can overlap.
The reference below is the former all-pairs loop; for random ORF sets,
including identical, touching, contained and wrap-sized intervals, the
records and their order must be identical.
"""

import random
import unittest
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from analyses._orf.summary import analyze_orf_relationships
from analyses._orf import codes
from core import api


def reference_relationships(orfs):
    relationships = []
    for i in range(len(orfs)):
        a = orfs[i]
        a_start, a_end = min(a["start"], a["stop"]), max(a["start"], a["stop"])
        for j in range(i + 1, len(orfs)):
            b = orfs[j]
            b_start, b_end = min(b["start"], b["stop"]), max(b["start"], b["stop"])
            overlap_start, overlap_end = max(a_start, b_start), min(a_end, b_end)
            if overlap_start > overlap_end:
                continue
            length = overlap_end - overlap_start + 1
            a_contains_b = a_start <= b_start and b_end <= a_end
            b_contains_a = b_start <= a_start and a_end <= b_end
            nesting = None
            if a_contains_b and not b_contains_a:
                relation, nesting = "Nested", (a, b)
            elif b_contains_a and not a_contains_b:
                relation, nesting = "Nested", (b, a)
            else:
                relation = "Overlapping"
            record = {"orf_a_rank": a["rank"], "orf_b_rank": b["rank"], "relation": relation,
                      "overlap_start": overlap_start + 1, "overlap_end": overlap_end + 1,
                      "overlap_length": length, "overlap_fraction_a": length / a["length"],
                      "overlap_fraction_b": length / b["length"],
                      "same_strand": a["strand"] == b["strand"], "same_frame": a["frame"] == b["frame"]}
            if nesting is not None:
                record["outer_orf_rank"] = nesting[0]["rank"]
                record["inner_orf_rank"] = nesting[1]["rank"]
            relationships.append(record)
    return relationships


def random_orfs(rng, count, span):
    orfs = []
    for rank in range(1, count + 1):
        start = rng.randrange(span)
        length = rng.choice([3, 30, 90, 300, 1500, span // 2])
        stop = min(span - 1, start + length - 1)
        strand = rng.choice("+-")
        if strand == "-" and rng.random() < 0.5:
            start, stop = stop, start          # either coordinate order occurs
        orfs.append({"rank": rank, "start": start, "stop": stop, "length": abs(stop - start) + 1,
                     "strand": strand, "frame": rng.choice([1, 2, 3])})
    if orfs:
        orfs.append(dict(orfs[0], rank=count + 1))   # an identical interval
    return orfs


class RelationshipSweepTests(unittest.TestCase):
    def test_random_sets_match_the_all_pairs_loop(self):
        rng = random.Random(20261120)
        for count, span in ((0, 10), (1, 10), (5, 50), (60, 500), (300, 20_000), (400, 2_000)):
            orfs = random_orfs(rng, count, span)
            with self.subTest(count=count, span=span):
                result = analyze_orf_relationships(orfs)
                expected = reference_relationships(orfs)
                self.assertEqual(result["relationships"], expected)
                self.assertEqual(result["relationship_count"], len(expected))
                self.assertEqual(result["nested_count"], sum(r["relation"] == "Nested" for r in expected))

    def test_real_orf_search_matches_the_all_pairs_loop(self):
        fasta = Path(__file__).parent / "fixtures" / "pUC19_L09137.2.fasta"
        sequence = "".join(line.strip() for line in fasta.read_text().splitlines() if not line.startswith(">"))
        result = api.find_orfs(sequence, min_length=30, start_codon_policy="atg_and_alternative",
                               nested=True, top_n=10_000)
        orfs = result.data["orfs"]
        self.assertGreater(len(orfs), 20)
        self.assertEqual(analyze_orf_relationships(orfs)["relationships"], reference_relationships(orfs))


class CodonExpansionCacheTests(unittest.TestCase):
    def test_callers_get_a_fresh_list(self):
        first = codes.expand_iupac_codon("ATN")
        first.append("XXX")
        self.assertEqual(codes.expand_iupac_codon("ATN"), ["ATA", "ATT", "ATC", "ATG"])
        self.assertEqual(codes.expand_iupac_codon("atg"), ["ATG"])
        self.assertEqual(codes.expand_iupac_codon("AT"), [])
        self.assertEqual(codes.expand_iupac_codon("AT?"), [])


if __name__ == "__main__":
    unittest.main()
