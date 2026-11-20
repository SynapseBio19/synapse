"""HGVS notation benchmarked against ClinVar (curated reference data).

The fixture tests/fixtures/clinvar/hgvs_benchmark.json holds ClinVar
variant names (e.g. ``NM_000546.6(TP53):c.215C>G (p.Pro72Arg)``) for six
genes and the RefSeq coding sequences they refer to; it was built by
scripts/benchmarks/clinvar_hgvs.py (source, date and selection are in the
file). ClinVar and RefSeq are public-domain NCBI records.

For every variant the observed coding sequence is rebuilt from ClinVar's
c. name. Synapse then compares reference and observed sequence as a user
would (global alignment, mutation calling, HGVS formatting with the 3'
rule), and its c. description must equal ClinVar's. For single-codon
substitutions the predicted protein change must equal ClinVar's p. name
after converting three-letter to one-letter codes (Synapse writes the
predicted form in parentheses, e.g. ``p.(P72R)``).

Result on the full set (1,500 variants, 2026-11-20): 1,489 c. names are
identical to ClinVar's. The other 11 are deletion-insertions that the
global alignment (linear gap costs) splits into several changes; Synapse
then writes an allele such as ``c.[217_218del;220_221del]`` for ClinVar's
``c.217_221delinsT``. Both describe the same observed sequence, which is
checked, but ClinVar's single delins is the preferred HGVS form. All 1,186
protein changes of codon substitutions agree.

By default the test evaluates every 10th variant and the 11 documented
cases (about 15 s); ``SYNAPSE_CLINVAR_FULL=1`` evaluates all 1,500
(about 2.5 min).

This checks notation and coordinate handling on real, curated variants.
It says nothing about pathogenicity, which Synapse does not assess.
"""

import hashlib
import json
import os
import re
import unittest
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import mutation as M
from tests.test_hgvs_allele_shift import apply_allele

FIXTURE = Path(__file__).parent / "fixtures" / "clinvar" / "hgvs_benchmark.json"
FIXTURE_SHA256 = "0348314292321b9646ff921505625a92b461e724718f53b18578d74e8df4a0c4"
FULL = os.environ.get("SYNAPSE_CLINVAR_FULL") == "1"

# ClinVar names whose single delins the alignment splits into several
# changes (see the module docstring); Synapse's allele is equivalent.
SPLIT_DELINS = {
    "c.221_227delinsCTCGGTGCCA", "c.217_221delinsT", "c.385_388delinsCCACA",
    "c.838_848delinsCA", "c.271_272delinsTTT", "c.975_988delinsCGCTT",
    "c.528_534delinsA", "c.1165_1169delinsGAGT", "c.1194_1222delinsCCA",
    "c.1197_1241delinsA", "c.1301_1325delinsAGCGGCCG",
}

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
    "Sec": "U",
}
C_NAME = re.compile(r"^c\.(\d+)(?:_(\d+))?(?:([ACGT])>([ACGT])|(del)|(dup)|ins([ACGT]+)|delins([ACGT]+))$")
P_SUBSTITUTION = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\?)$")


def apply_c_name(cds, name):
    """Observed coding sequence for a ClinVar c. name (exonic changes only)."""
    match = C_NAME.match(name)
    if match is None:
        raise ValueError(name)
    first = int(match.group(1))
    last = int(match.group(2) or first)
    ref_base, alt_base, deletion, duplication, insertion, delins = match.groups()[2:]
    if ref_base:
        if cds[first - 1] != ref_base:
            raise ValueError(f"{name}: reference base is {cds[first - 1]}")
        return cds[:first - 1] + alt_base + cds[first:]
    if deletion:
        return cds[:first - 1] + cds[last:]
    if duplication:
        return cds[:last] + cds[first - 1:last] + cds[last:]
    if insertion:
        if last != first + 1:
            raise ValueError(name)
        return cds[:first] + insertion + cds[first:]
    return cds[:first - 1] + delins + cds[last:]


def clinvar_protein_to_synapse(p_name):
    """ClinVar's p. name in Synapse's form: p.Pro72Arg -> p.(P72R).

    Two ClinVar shorthands are mapped to the HGVS recommendation that
    Synapse follows (den Dunnen et al. 2016): a changed initiation codon
    (ClinVar p.Met1Arg) is p.(M1?), because the effect on initiation is
    unknown, and a lost stop codon (ClinVar p.Ter214Leu) is an extension
    with an undetermined new stop, p.(*214Lext*?). None for forms that are
    not a single-residue change.
    """
    if p_name is None:
        return None
    match = P_SUBSTITUTION.match(p_name)
    if match is None or match.group(1) not in THREE_TO_ONE:
        return None
    ref, position, alt = match.groups()
    alt = alt if alt in ("=", "?") else THREE_TO_ONE.get(alt)
    if alt is None:
        return None
    ref = THREE_TO_ONE[ref]
    if ref == "M" and position == "1" and alt != "=":
        return "p.(M1?)"
    if ref == "*" and alt not in ("*", "="):
        return f"p.(*{position}{alt}ext*?)"
    return f"p.({ref}{position}{alt})"


def is_in_frame_substitution(c_name):
    """A c. change that replaces bases one for one (X>Y, or delins of equal length)."""
    match = C_NAME.match(c_name)
    assert match is not None, c_name
    if match.group(3):
        return True
    if match.group(8):
        first, last = int(match.group(1)), int(match.group(2) or match.group(1))
        return len(match.group(8)) == last - first + 1
    return False


def synapse_descriptions(cds, observed):
    result = M.compare_sequences(cds, observed)
    mutations = list(result.mutations)
    c_name = M.hgvs_c_allele(mutations, cds)
    # Changes inside one codon are called as several mutations that carry
    # the same codon consequence. A change across two codons where one codon
    # stays synonymous (KRAS c.36_37delinsGC: p.(G12=) and p.(G13R)) is one
    # protein change at the protein level, as ClinVar writes it.
    p_names = {name for mutation in mutations if (name := M.mutation_hgvs_p(mutation)) is not None}
    if len(p_names) > 1:
        p_names = {name for name in p_names if not name.endswith("=)")}
    p_name = p_names.pop() if len(p_names) == 1 else None
    return c_name, p_name


def selected_variants(fixture, full=FULL):
    variants = fixture["variants"]
    if full:
        return list(variants)
    return [v for index, v in enumerate(variants) if index % 10 == 0 or v["c"] in SPLIT_DELINS]


def evaluate(fixture, full=FULL):
    """Return per-variant outcomes (shared by the tests and docs/BENCHMARKS.md)."""
    outcomes = []
    for variant in selected_variants(fixture, full):
        cds = fixture["transcripts"][variant["transcript"]]["cds"]
        observed = apply_c_name(cds, variant["c"])
        c_name, p_name = synapse_descriptions(cds, observed)
        # Synapse writes a protein change for codon substitutions only, not
        # frameshift or in-frame indel notation (docs/METHODS.md).
        expected_p = (clinvar_protein_to_synapse(variant["p"])
                      if is_in_frame_substitution(variant["c"]) else None)
        outcomes.append({
            "uid": variant["uid"], "gene": variant["gene"], "clinvar_c": variant["c"],
            "synapse_c": c_name, "clinvar_p": variant["p"], "expected_p": expected_p,
            "synapse_p": p_name,
            "c_equal": c_name == variant["c"],
            "rebuilds_observed": apply_allele(cds, c_name) == observed,
            "p_checked": expected_p is not None,
            "p_equal": expected_p is None or p_name == expected_p,
        })
    return outcomes


@unittest.skipUnless(FIXTURE.exists(), "ClinVar fixture not built (scripts/benchmarks/clinvar_hgvs.py)")
class ClinVarHgvsBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
        cls.outcomes = evaluate(cls.fixture)

    def test_fixture_is_the_documented_set(self):
        self.assertEqual(hashlib.sha256(FIXTURE.read_bytes()).hexdigest(), FIXTURE_SHA256)
        self.assertEqual(set(self.fixture["transcripts"]),
                         {"NM_000518.5", "NM_004985.5", "NM_000551.4", "NM_000546.6", "NM_000314.8", "NM_001110792.2"})
        for accession, entry in self.fixture["transcripts"].items():
            with self.subTest(transcript=accession):
                self.assertEqual(len(entry["cds"]), entry["cds_length"])
                self.assertEqual(entry["cds"][:3], "ATG")
                self.assertEqual(len(entry["cds"]) % 3, 0)
        self.assertEqual(len(self.fixture["variants"]), 1500)
        self.assertEqual(SPLIT_DELINS - {v["c"] for v in self.fixture["variants"]}, set())

    def test_c_names_match_clinvar_except_the_split_delins(self):
        different = {o["clinvar_c"] for o in self.outcomes if not o["c_equal"]}
        self.assertEqual(different, SPLIT_DELINS)

    def test_every_description_rebuilds_the_observed_sequence(self):
        self.assertEqual([o["clinvar_c"] for o in self.outcomes if not o["rebuilds_observed"]], [])

    def test_every_protein_substitution_matches_clinvar(self):
        checked = [o for o in self.outcomes if o["p_checked"]]
        self.assertGreaterEqual(len(checked), 1186 if FULL else 110)
        different = [(o["gene"], o["clinvar_c"], o["clinvar_p"], o["synapse_p"]) for o in checked if not o["p_equal"]]
        self.assertEqual(different, [])

    def test_all_change_types_are_represented(self):
        kinds = {"substitution": r">", "deletion": r"\d(del)$", "duplication": r"dup$",
                 "insertion": r"\d(ins)[ACGT]+$", "deletion-insertion": r"delins"}
        for kind, pattern in kinds.items():
            with self.subTest(kind=kind):
                self.assertTrue(any(re.search(pattern, o["clinvar_c"]) for o in self.outcomes))


if __name__ == "__main__":
    unittest.main()
