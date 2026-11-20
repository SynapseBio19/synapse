"""Author names merged from several sources are listed once (2026-11-20).

PubMed gives full names ("Ann Smith"); Europe PMC's authorString gives
the family name followed by initials without a comma ("Smith A"). The
de-duplication read "Smith A" as given name "Smith" and family name "A",
so a merged record listed the same author twice, in the result and in
every export (found with the literature characterization scenarios).

A trailing token of one to three capital initials is now also read as
"family name, initials". The first reading is kept as well, so names
whose family name is written in capitals ("Wei LI") still match their
other spellings.
"""

import unittest

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import literature as L


class AuthorMergeTests(unittest.TestCase):
    def check(self, names, expected):
        self.assertEqual(L._deduplicate_author_names(names), expected)

    def test_family_name_with_initials_matches_the_full_name(self):
        self.check(["Ann Smith", "Smith A"], ["Ann Smith"])
        self.check(["Smith A", "Ann Smith"], ["Ann Smith"])           # the longer form is kept
        self.check(["Anna Beth Smith", "Smith AB"], ["Anna Beth Smith"])
        self.check(["Test Consortium", "Consortium T."], ["Test Consortium"])

    def test_capitalized_family_names_still_match(self):
        self.check(["Wei LI", "Wei Li"], ["Wei LI"])
        self.check(["Wei Li", "Li W"], ["Wei Li"])

    def test_distinct_people_stay_distinct(self):
        self.check(["John Smith", "Jane Smith"], ["John Smith", "Jane Smith"])
        self.check(["Ann Smith", "Smith B"], ["Ann Smith", "Smith B"])
        self.check(["Ann Smith", "Jones A"], ["Ann Smith", "Jones A"])

    def test_former_behavior_is_unchanged(self):
        # The built-in self-test's cases.
        self.check(["A. Krull", "Abigail G Krull", "Seth Entriken", "Seth M Entriken"],
                   ["Abigail G Krull", "Seth Entriken", "Seth M Entriken"])
        self.check(["A. Krull", "Abigail G Krull", "Seth M Entriken", "Seth Entriken"],
                   ["Abigail G Krull", "Seth M Entriken", "Seth Entriken"])
        self.check(["Ann Smith", "A. Smith", "Bo Li", "bo li"], ["Ann Smith", "Bo Li"])

    def test_merged_record_lists_the_author_once(self):
        pubmed = L.LiteratureArticle(title="Lactase persistence in Europe", doi="10.1000/abc",
                                     authors=["Ann Smith", "Test Consortium"])
        europe = L.LiteratureArticle(title="Lactase persistence in Europe.", doi="10.1000/ABC",
                                     authors=["Smith A", "Consortium T."])
        merged = L.deduplicate_articles([pubmed, europe])
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0].authors, ["Ann Smith", "Test Consortium"])


if __name__ == "__main__":
    unittest.main()
