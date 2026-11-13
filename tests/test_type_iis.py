"""P2-7: Type IIS enzymes cut outside an asymmetric site, on either strand.

Worked by hand (BsaI, GGTCTC(1/5)):

    index   0123 456789 0123456789 0123 456789 0123
    top     AAAA GGTCTC AAAAAAAAAA TTTT GAGACC CCCC

* "+" site at 4: top cut 4 + 7 = 11, bottom cut 4 + 11 = 15, 5' overhang
  top[11:15] = "AAAA".
* "-" site (GAGACC = reverse complement of GGTCTC) at 24: the cuts lie on
  the other side, top 24 + 6 - 11 = 19 and bottom 24 + 6 - 7 = 23, overhang
  top[19:23] = "ATTT".
"""

from __future__ import annotations

import unittest

from tests import _support  # noqa: F401
from ui import restriction_menu, style

SEQUENCE = "AAAA" + "GGTCTC" + "A" * 10 + "TTTT" + "GAGACC" + "CCCC"


class TypeIISSiteTests(unittest.TestCase):
    def test_both_orientations_with_their_cuts_and_overhangs(self):
        from analyses.restriction import find_type_iis_sites
        sites = find_type_iis_sites(SEQUENCE, "BsaI")
        self.assertEqual([(s["strand"], s["site_position_0based"], s["top_cut_0based"], s["bottom_cut_0based"],
                           s["overhang"]) for s in sites],
                         [("+", 4, 11, 15, "AAAA"), ("-", 24, 19, 23, "ATTT")])
        # unsuffixed keys are 1-based, as in the restriction map: site start, first base after the cut
        self.assertEqual([(s["site_position"], s["top_cut"], s["bottom_cut"]) for s in sites],
                         [(5, 12, 16), (25, 20, 24)])
        self.assertTrue(all(s["status"] == "CONFIRMED" and s["cuts_inside"] for s in sites))

    def test_digest_uses_both_strands(self):
        from analyses.restriction import perform_digest
        digest = perform_digest(SEQUENCE, ["BsaI"], "linear")
        self.assertEqual(digest["cut_positions"], [11, 19])
        self.assertEqual([f["length"] for f in digest["fragments"]], [11, 8, 15])
        self.assertEqual({c["strand"] for c in digest["cuts"]}, {"+", "-"})

    def test_a_cut_beyond_a_linear_end_does_not_count(self):
        from analyses.restriction import find_type_iis_sites, perform_digest
        near_end = "ACGTACGTAC" + "GGTCTC" + "AAA"      # top cut would fall at 17 of 19, bottom at 21
        site = find_type_iis_sites(near_end, "BsaI")[0]
        self.assertFalse(site["cuts_inside"])
        self.assertEqual(perform_digest(near_end, ["BsaI"], "linear")["cut_count"], 0)

    def test_circular_cuts_wrap_around_the_origin(self):
        from analyses.restriction import find_type_iis_sites, perform_digest
        sequence = "TTAACCGG" + "A" * 30 + "GGTCTC"   # 44 bp; cuts at 45 and 49 wrap to 1 and 5
        site = find_type_iis_sites(sequence, "BsaI", topology="circular")[0]
        self.assertTrue(site["cuts_inside"])
        self.assertEqual((site["top_cut_0based"], site["bottom_cut_0based"]), (1, 5))
        self.assertEqual(site["overhang"], "TAAC")
        self.assertTrue(site["wraps_origin"])
        self.assertEqual(perform_digest(sequence, ["BsaI"], "circular")["cut_positions"], [1])
        self.assertEqual(perform_digest(sequence, ["BsaI"], "linear")["cut_count"], 0)

    def test_an_ambiguous_site_is_possible_not_a_cut(self):
        from analyses.restriction import find_type_iis_sites, perform_digest
        sequence = "A" * 10 + "GGTNTC" + "A" * 20
        self.assertEqual([s["status"] for s in find_type_iis_sites(sequence, "BsaI")], ["POSSIBLE"])
        digest = perform_digest(sequence, ["BsaI"], "linear")
        self.assertEqual((digest["cut_count"], len(digest["possible_cuts"])), (0, 1))

    def test_three_nucleotide_overhang_of_sapi(self):
        from analyses.restriction import find_type_iis_sites
        sequence = "A" * 5 + "GCTCTTC" + "TCAGGG" + "A" * 10
        site = find_type_iis_sites(sequence, "SapI")[0]
        self.assertEqual((site["top_cut_0based"], site["bottom_cut_0based"], site["overhang"]), (13, 16, "CAG"))

    def test_type_iis_enzymes_stay_out_of_the_fixed_overhang_cloning_model(self):
        from analyses import cloning
        from data.restriction_enzymes import RESTRICTION_ENZYMES, TYPE_IIS_ENZYMES
        self.assertFalse(set(TYPE_IIS_ENZYMES) & set(RESTRICTION_ENZYMES))
        with self.assertRaises(Exception):
            cloning.get_overhang("BsaI")

    def test_digest_dialog_lists_and_reports_type_iis_cuts(self):
        import builtins
        import contextlib
        import io
        from unittest import mock
        answers = iter(["BsaI", "", ""])
        buffer = io.StringIO()
        with mock.patch.object(builtins, "input", lambda prompt="": next(answers)),                 contextlib.redirect_stdout(buffer):
            result = restriction_menu._restriction_digest_gel_ui(SEQUENCE)
        lines = buffer.getvalue().splitlines()
        self.assertIn("23. BsaI       GGTCTC(1/5)", lines)
        self.assertIn("TYPE IIS CUTS", lines)
        widths = (8, 7, 10, 10, 13)
        self.assertIn(style.columns(("BsaI", "+", 5, "after 11", "after 15", "AAAA"), widths), lines)
        self.assertIn(style.columns(("BsaI", "-", 25, "after 19", "after 23", "ATTT"), widths), lines)
        self.assertEqual(result["digest"]["cut_positions"], [11, 19])

    def test_digest_events_follow_the_map_convention(self):
        """site_position and top/bottom cuts are 1-based in maps and digests alike (2026-11-11)."""
        from analyses.restriction import analyze_restriction_site, perform_digest
        sequence = "A" * 20 + "GAATTC" + "A" * 20
        (event,) = perform_digest(sequence, ["EcoRI"], "linear")["cuts"]
        mapped = analyze_restriction_site(sequence, "EcoRI", 20)
        for key in ("site_position", "site_position_0based", "top_cut", "top_cut_0based",
                    "bottom_cut", "bottom_cut_0based"):
            with self.subTest(key=key):
                self.assertEqual(event[key], mapped[key])
        self.assertEqual((event["site_position"], event["top_cut"], event["position_0based"]), (21, 22, 21))

    def test_unknown_names_are_still_rejected(self):
        from analyses.restriction import find_type_iis_sites, perform_digest
        with self.assertRaises(ValueError):
            perform_digest(SEQUENCE, ["BsaX"], "linear")
        with self.assertRaises(ValueError):
            find_type_iis_sites(SEQUENCE, "EcoRI")


if __name__ == "__main__":
    unittest.main()
