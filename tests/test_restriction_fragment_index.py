"""Faster restriction digests give exactly the former results (2026-11-20).

1. Fragment annotation. 
perform_digest stores, for each fragment, the cut events at its two
boundaries. Until 2026-11-20 every fragment scanned the whole cut list
(fragments x cuts); with all 22 enzymes on the 4.6 Mb E. coli genome the
digest did not finish within 15 minutes (docs/BENCHMARKS.md). The cuts are
now indexed by coordinate. The reference below is the former loop
(linear, and circular with its explicit wrap clause); on real and random
sequences, both topologies and several enzyme panels the digest must be
identical.

2. Recognition-site search. Windows were tested base by base with
bases_match; a regular expression with the same compatibility classes now
finds them. The reference is the former per-position loop, on random
sequences with IUPAC ambiguity codes and lower-case letters, linear and
circular, for every database recognition sequence and random IUPAC ones.
"""

import random
import unittest
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from analyses._restriction.digest import perform_digest
from analyses._restriction.sites import IUPAC_BASES, _restriction_sequence_status, _scan_recognition_sites, bases_match
from data.restriction_enzymes import RESTRICTION_ENZYMES, TYPE_IIS_ENZYMES

FIXTURES = Path(__file__).parent / "fixtures"


def reference_annotation(digest):
    """The fragments' cut_events/enzymes recomputed with the former loops."""
    cuts = digest["cuts"]
    annotated = []
    for fragment in digest["fragments"]:
        boundaries = {fragment["start_0based"], fragment["end_0based_exclusive"]}
        if digest["topology"] == "linear":
            events = [cut for cut in cuts if cut["position_0based"] in boundaries]
        else:
            events = [cut for cut in cuts
                      if cut["position_0based"] in boundaries
                      or (fragment["wraps_origin"] and cut["position_0based"] == 0 and 0 in boundaries)]
        annotated.append((
            [{"enzyme": cut["enzyme"], "position_1based": cut["position_1based"]} for cut in events],
            sorted({cut["enzyme"] for cut in events}),
        ))
    return annotated


class FragmentIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        lines = (FIXTURES / "pUC19_L09137.2.fasta").read_text().splitlines()
        puc19 = "".join(line.strip() for line in lines if line and not line.startswith(">"))
        rng = random.Random(20261120)
        cls.sequences = [puc19, "".join(rng.choice("ACGT") for _ in range(6000)),
                         "GAATTC" + "".join(rng.choice("ACGT") for _ in range(400)) + "GAATTC"]

    def test_annotation_matches_the_former_loops(self):
        panels = [sorted(RESTRICTION_ENZYMES), ["EcoRI"], ["AluI", "HaeIII", "MspI", "HpaII"],
                  sorted(TYPE_IIS_ENZYMES) + ["EcoRI", "TaqI"]]
        checked = 0
        for sequence in self.sequences:
            for topology in ("linear", "circular"):
                for panel in panels:
                    with self.subTest(length=len(sequence), topology=topology, panel=panel[:3]):
                        digest = perform_digest(sequence, panel, topology=topology)
                        got = [(f["cut_events"], f["enzymes"]) for f in digest["fragments"]]
                        self.assertEqual(got, reference_annotation(digest))
                        checked += sum(1 for events, _ in got if events)
        self.assertGreater(checked, 500)

    def test_shared_coordinates_keep_the_cut_list_order(self):
        # MspI and HpaII cut CCGG at the same coordinate: both are events.
        digest = perform_digest("AAAACCGGAAAA", ["MspI", "HpaII"], topology="linear")
        self.assertEqual(digest["fragments"][0]["enzymes"], ["HpaII", "MspI"])
        self.assertEqual([e["enzyme"] for e in digest["fragments"][0]["cut_events"]],
                         [cut["enzyme"] for cut in digest["cuts"]])


def reference_scan(sequence, recognition, circular):
    """The former _scan_recognition_sites loop."""
    matches = []
    search_space, last_start = sequence, len(sequence) - len(recognition)
    if circular and len(sequence) >= len(recognition):
        search_space = sequence + sequence[:len(recognition) - 1]
        last_start = len(sequence) - 1
    for position in range(last_start + 1):
        fragment = search_space[position:position + len(recognition)]
        if not all(bases_match(a, b) for a, b in zip(fragment, recognition)):
            continue
        offsets = [i + 1 for i, base in enumerate(fragment) if base not in {"A", "C", "G", "T"}]
        matches.append({"position": position, "position_1based": position + 1,
                        "status": "POSSIBLE" if offsets else "CONFIRMED", "ambiguous_offsets": offsets,
                        "sequence": fragment, "recognition_sequence": recognition,
                        "wraps_origin": position + len(recognition) > len(sequence)})
    return tuple(matches)


class SiteScanTests(unittest.TestCase):
    def test_regex_scan_matches_the_base_by_base_loop(self):
        rng = random.Random(11)
        recognitions = sorted({e["recognition_sequence"] for e in RESTRICTION_ENZYMES.values()}
                              | {e["recognition_sequence"] for e in TYPE_IIS_ENZYMES.values()})
        recognitions += ["".join(rng.choice("ACGTRYSWKMBDHVN") for _ in range(rng.randint(1, 8)))
                         for _ in range(25)]
        alphabet = "ACGT" * 12 + "RYSWKMBDHVN" + "acgtn"
        sequences = ["".join(rng.choice(alphabet) for _ in range(n)) for n in (1, 5, 300, 3000)]
        sequences.append("GAATTCNNGAATTCGGATCC")
        found = 0
        for sequence in sequences:
            for recognition in recognitions:
                for circular in (False, True):
                    with self.subTest(length=len(sequence), site=recognition, circular=circular):
                        expected = reference_scan(sequence, recognition, circular)
                        self.assertEqual(_scan_recognition_sites(sequence, recognition, circular), expected)
                        found += len(expected)
        self.assertGreater(found, 1000)
        self.assertEqual(_scan_recognition_sites("ACGT", "AXG", False), reference_scan("ACGT", "AXG", False))


def reference_status(sequence):
    """The former _restriction_sequence_status (position lists always enumerated)."""
    normalized = "".join(sequence.upper().split()) if isinstance(sequence, str) else ""
    invalid = [i + 1 for i, base in enumerate(normalized) if base not in IUPAC_BASES]
    return {"sequence": normalized, "valid": bool(normalized) and not invalid, "length": len(normalized),
            "invalid_positions": invalid, "invalid_symbols": sorted({normalized[p - 1] for p in invalid}),
            "ambiguous_positions": [i + 1 for i, base in enumerate(normalized)
                                    if base not in {"A", "C", "G", "T"} and base in IUPAC_BASES]}


class SequenceStatusTests(unittest.TestCase):
    def test_status_shortcut_matches_the_former_enumeration(self):
        rng = random.Random(5)
        for alphabet in ("ACGT", "ACGTN", "ACGTX-", "acgtryn", "AC GT	"):
            for length in (0, 1, 50, 500):
                sequence = "".join(rng.choice(alphabet) for _ in range(length))
                with self.subTest(alphabet=alphabet, length=length):
                    self.assertEqual(_restriction_sequence_status(sequence), reference_status(sequence))
        self.assertEqual(_restriction_sequence_status(None), reference_status(None))


if __name__ == "__main__":
    unittest.main()
