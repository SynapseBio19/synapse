"""P1-9: ORFs across the origin of a circular sequence and 5'-partial ORFs.

Circular: a plasmid has no ends, so an ORF may start before the arbitrary
origin and stop after it. The search reads the circle as dna + dna[:n-1];
coordinates stay unrolled (start <= stop) and stop_on_circle gives the last
base on the circle. On a circle the set of ORFs must not depend on where the
origin lies, which the rotation test checks directly.

5'-partial: in a linear fragment the start codon of an ORF may lie upstream
of the sequence. Such an open frame runs in from the 5' end and is reported,
on request, as "Partial - Missing Start"; its first codon is not a start
codon and is not read as initiator methionine.
"""

from __future__ import annotations

import random
import unittest

from tests import _support  # noqa: F401
from tests._support import reverse_complement
from ui import orf_menu

BODY = "AAA" * 30
# 6 + 3 + 20 + 3 + 90 = 122 bp: ATG at 29, 30 AAA codons to the end, then across
# the origin AAA AAA and the stop TAA at 6-8.
CIRCLE = BODY[:6] + "TAA" + "C" * 20 + "ATG" + BODY


class OriginSpanningOrfTests(unittest.TestCase):
    def test_an_orf_runs_across_the_origin(self):
        from analyses.orf import find_orfs, orf_end_on_sequence
        (orf,) = [o for o in find_orfs(CIRCLE, min_length=30, topology="circular") if o["strand"] == "+"]
        self.assertEqual((orf["start"], orf["stop"], orf["stop_on_circle"]), (29, 130, 8))
        self.assertEqual((orf["length"], orf["status"], orf["wraps_origin"]), (102, "Complete", True))
        self.assertEqual(orf["sequence"], "ATG" + BODY + "AAAAAA" + "TAA")
        self.assertEqual(orf_end_on_sequence(orf), 9)

    def test_linear_reading_keeps_the_partial_and_no_wrap(self):
        from analyses.orf import find_orfs
        plus = [o for o in find_orfs(CIRCLE, min_length=30) if o["strand"] == "+"]
        self.assertEqual([(o["start"], o["status"]) for o in plus], [(29, "Partial - Missing Stop")])
        self.assertFalse(any(o.get("wraps_origin") for o in plus))

    def test_the_minus_strand_wraps_too(self):
        from analyses.orf import find_orfs
        minus = [o for o in find_orfs(reverse_complement(CIRCLE), min_length=30, topology="circular")
                 if o["strand"] == "-"]
        self.assertEqual(len(minus), 1)
        self.assertTrue(minus[0]["wraps_origin"])
        self.assertEqual(minus[0]["sequence"], "ATG" + BODY + "AAAAAA" + "TAA")

    def test_orfs_do_not_depend_on_where_the_origin_lies(self):
        from analyses.orf import find_orfs
        rng = random.Random(1909)
        circle = "".join(rng.choice("ACGT") for _ in range(900))

        def found(sequence):
            return sorted((o["strand"], o["sequence"]) for o in
                          find_orfs(sequence, min_length=60, topology="circular"))

        reference = found(circle)
        self.assertTrue(any(True for _ in reference))
        for rotation in (1, 2, 450, 899):
            with self.subTest(rotation=rotation):
                self.assertEqual(found(circle[rotation:] + circle[:rotation]), reference)

    def test_the_longest_orf_per_stop_wins_across_the_origin(self):
        from analyses.orf import find_orfs
        # a second ATG after the origin, in frame with the wrapping ORF and before its stop
        circle = "ATGAAA" + "TAA" + "C" * 20 + "ATG" + BODY
        plus = [o for o in find_orfs(circle, min_length=9, topology="circular") if o["strand"] == "+"]
        ends = [o["stop_on_circle"] if o.get("wraps_origin") else o["stop"] for o in plus]
        self.assertEqual(len(ends), len(set(ends)))            # one ORF per stop
        self.assertTrue(any(o.get("wraps_origin") for o in plus))
        nested = [o for o in find_orfs(circle, min_length=9, topology="circular", nested=True)
                  if o["strand"] == "+"]
        self.assertGreater(len(nested), len(plus))              # the inner ATG is kept when nested


class FivePrimePartialOrfTests(unittest.TestCase):
    # frame +1: AAA GGG CCC ... AAA TAG (stop at 30-32), no start codon before it
    SEQUENCE = "AAAGGGCCCAAAGGGCCCAAAGGGCCCAAATAGCCC"

    def test_off_by_default(self):
        from analyses.orf import find_orfs
        self.assertFalse(any(o["status"].startswith("Partial - Missing Start")
                             for o in find_orfs(self.SEQUENCE, min_length=30)))

    def test_an_open_frame_from_the_5_prime_end(self):
        from analyses.orf import find_orfs
        orfs = find_orfs(self.SEQUENCE, min_length=30, include_5prime_partial=True, include_partial=False)
        (orf,) = [o for o in orfs if o["frame"] == 1]
        self.assertEqual((orf["start"], orf["stop"], orf["length"]), (0, 32, 33))
        self.assertEqual((orf["status"], orf["start_status"], orf["stop_codon"]),
                         ("Partial - Missing Start", "absent", "TAG"))
        self.assertEqual(orf["protein_length"], 10)

    def test_a_frame_that_begins_with_a_start_codon_is_complete(self):
        from analyses.orf import find_orfs
        orfs = find_orfs("ATG" + "AAA" * 10 + "TAA", min_length=30, include_5prime_partial=True)
        self.assertEqual([o["status"] for o in orfs if o["frame"] == 1], ["Complete"])

    def test_complete_orfs_inside_the_open_frame_are_kept(self):
        from analyses.orf import find_orfs
        sequence = "GCC" * 10 + "ATG" + "AAA" * 10 + "TAA"
        frame_one = [o["status"] for o in find_orfs(sequence, min_length=30, include_5prime_partial=True)
                     if o["frame"] == 1]
        self.assertEqual(sorted(frame_one), ["Complete", "Partial - Missing Start"])

    def test_the_first_codon_is_not_read_as_initiator(self):
        from analyses.orf import find_orfs, translate_orf_sequence
        from core import api
        (orf,) = [o for o in find_orfs("TTGAAAAAAAAAAAATAA", min_length=12, include_5prime_partial=True)
                  if o["status"] == "Partial - Missing Start"]
        self.assertEqual(translate_orf_sequence(orf), "LKKKK")       # TTG is Leu here, not Met
        result = api.analyze_protein_from_cds(orf)
        self.assertEqual(result.data["protein"], "LKKKK")

    def test_minus_strand_coordinates(self):
        from analyses.orf import find_orfs
        minus = reverse_complement(self.SEQUENCE)
        (orf,) = [o for o in find_orfs(minus, min_length=30, include_5prime_partial=True, include_partial=False)
                  if o["strand"] == "-" and o["frame"] == -1]
        length = len(minus)
        self.assertEqual((orf["start"], orf["stop"]), (length - 33, length - 1))
        self.assertEqual(orf["sequence"], self.SEQUENCE[:33])


class OrfOptionsReachTheResultsTests(unittest.TestCase):
    def test_parameters_record_topology_and_partials(self):
        from analyses.orf import analyze_orfs
        params = analyze_orfs(CIRCLE, min_length=30, topology="circular")["parameters"]
        self.assertEqual((params["topology"], params["include_5prime_partial"]), ("circular", False))

    def test_pipeline_uses_the_record_topology(self):
        from core.pipeline import SynapsePipeline
        pipeline = SynapsePipeline(CIRCLE, "plasmid", topology="circular")
        result = pipeline.run_orf(min_length=30)
        self.assertTrue(any(o.get("wraps_origin") for o in result["orfs"]))

    def test_api_states_the_circular_reading(self):
        from core import api
        from core.models import SequenceRecord
        result = api.find_orfs(SequenceRecord("plasmid", CIRCLE, topology="circular"), min_length=30)
        self.assertTrue(any(o.get("wraps_origin") for o in result.data["orfs"]))
        self.assertTrue(any("across the origin" in text for text in result.limitations))

    def test_settings_cycle_through_the_partial_options(self):
        import builtins
        import contextlib
        import io
        from unittest import mock
        from analyses import orf as O
        answers = iter(["5", "5", "5", "0"])
        settings = {"min_length": 75, "start_codon_policy": "atg_only", "nested": False,
                    "ranking": "length", "include_partial": True, "genetic_code": 1}
        buffer = io.StringIO()
        with mock.patch.object(builtins, "input", lambda prompt="": next(answers)), \
                contextlib.redirect_stdout(buffer):
            result = orf_menu._orf_settings_menu(settings)
        shown = [line for line in buffer.getvalue().splitlines() if line.startswith("[5] ")]
        self.assertEqual([line.split(":", 1)[1].strip() for line in shown],
                         ["3' end (missing stop)", "5' and 3' ends", "excluded", "3' end (missing stop)"])
        self.assertEqual((result["include_partial"], result["include_5prime_partial"]), (True, False))


if __name__ == "__main__":
    unittest.main()
