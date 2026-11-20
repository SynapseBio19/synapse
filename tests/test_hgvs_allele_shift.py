"""The 3' rule must not move one change of an allele past the next one.

Found with the ClinVar benchmark (2026-11-20): for PTEN
NM_000314.8 c.528_534delinsA the alignment calls three deletions
(526_528, 531_532, 534). hgvs_c_allele shifted each of them to its most
3' position as if it were alone; the first one travelled through the
TATTATTAT repeat to 532_534 and the allele read
``c.[532_534del;531_532del;534del]``, overlapping deletions that describe
no sequence. HGVS applies the 3' rule to each variant without passing a
neighbouring variant, so a shift stops before the next change and the
changes of an allele stay separate. Every allele description must
rebuild the observed sequence.
"""

import json
import re
import unittest
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import mutation as M

FIXTURE = Path(__file__).parent / "fixtures" / "clinvar" / "hgvs_benchmark.json"
EVENT = re.compile(r"^(\d+)(?:_(\d+))?(?:([ACGT])>([ACGT])|(del)|(dup)|ins([ACGT]+)|delins([ACGT]+))$")


def events_of(allele):
    body = allele[2:]
    parts = body[1:-1].split(";") if body.startswith("[") else [body]
    events = []
    for part in parts:
        match = EVENT.match(part)
        if match is None:
            raise AssertionError(f"not an HGVS change: {part}")
        first = int(match.group(1))
        last = int(match.group(2) or first)
        events.append((first, last, match.groups()[2:]))
    return events


def apply_allele(reference, allele):
    """Apply a c.[...] allele; the changes must be in order and must not overlap."""
    events = events_of(allele)
    spans = []
    for first, last, kind in events:
        insertion = kind[4] is not None
        spans.append((first, last, insertion))
    for (a_first, a_last, a_ins), (b_first, b_last, b_ins) in zip(spans, spans[1:]):
        # An insertion is between its two flanking bases; others occupy first..last.
        a_end = a_first if a_ins else a_last
        b_start = b_last if b_ins else b_first
        if a_end >= b_start:
            raise AssertionError(f"overlapping or unordered changes in {allele}")
    sequence = reference
    for first, last, (ref, alt, deletion, duplication, insertion, delins) in reversed(events):
        if ref:
            sequence = sequence[:first - 1] + alt + sequence[first:]
        elif deletion:
            sequence = sequence[:first - 1] + sequence[last:]
        elif duplication:
            sequence = sequence[:last] + sequence[first - 1:last] + sequence[last:]
        elif insertion:
            sequence = sequence[:first] + insertion + sequence[first:]
        else:
            sequence = sequence[:first - 1] + delins + sequence[last:]
    return sequence


def allele(reference, observed):
    result = M.compare_sequences(reference, observed)
    return M.hgvs_c_allele(list(result.mutations), reference)


class AlleleShiftTests(unittest.TestCase):
    def test_synthetic_repeat_with_two_deletions(self):
        # Two deletions in one repeat: alone, the first would shift onto the second.
        reference = "ATGCCC" + "TATTATTAT" + "GCAGGCTAA"
        observed = "ATGCCC" + "TAT" + "TTAT" + "GCAGGCTAA"[1:]
        description = allele(reference, observed)
        self.assertEqual(apply_allele(reference, description), observed)

    def test_pten_case_from_clinvar(self):
        fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
        cds = fixture["transcripts"]["NM_000314.8"]["cds"]
        observed = cds[:527] + "A" + cds[534:]          # c.528_534delinsA
        description = allele(cds, observed)
        self.assertEqual(apply_allele(cds, description), observed)
        self.assertNotEqual(description, "c.[532_534del;531_532del;534del]")

    def test_a_single_change_still_takes_its_most_3_prime_position(self):
        reference = "ATGAAACGTTGGCTGCAAGGTTAA"
        self.assertEqual(allele(reference, "ATGAACGTTGGCTGCAAGGTTAA"), "c.6del")
        self.assertEqual(allele(reference, "ATGAAAACGTTGGCTGCAAGGTTAA"), "c.6dup")


if __name__ == "__main__":
    unittest.main()
