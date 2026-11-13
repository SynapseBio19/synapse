"""Variant normalization compared with the published output of `bcftools norm`.

The inputs, references and expected outputs are the bcftools 1.25 test
vectors (tests/fixtures/bcftools_norm/SOURCE.md). bcftools itself is not
run here; its own test suite asserts these outputs.

Synapse always splits multi-allelic records, so the split-then-normalize
case (`bcftools norm -m-`) is the like-for-like comparison. Joint
normalization of multi-allelic records (`norm.out`) is used only for the
REF check (`-cx`).
"""

import hashlib
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from tests import _support  # noqa: F401  (blocks network and input())
from analyses import variant_analysis as V

FIXTURES = Path(__file__).parent / "fixtures" / "bcftools_norm"
SHA256 = {
    "norm.fa": "da20982fbeb50362b580d8d2e06c492bfe221cb92020729a3e3071147764b060",
    "norm.out": "dba0b08005378b011f90695ee1ac5786eaf9b0caab767c6f7e28a54fd87a6b0a",
    "norm.right-align.1.out": "f16444cd7521694a32e8e18ce041afc63301800e200f7892d8def419c7722bc6",
    "norm.right-align.fa": "0eb5e4aefeabacac15326d354dbf4a478415a92f31d926332457e8048547d497",
    "norm.right-align.vcf": "85004daeb31afd43dba635150c8cea707e4963a21ae5f194a90e3018ff015a59",
    "norm.split.and.norm.out": "0e4412531059a5e49757a27a9d96e00dc33375301dcb4eb71b1f34e11110b221",
    "norm.split.vcf": "0363806f590a1ce92c7e2a8c3e16d976577e91fec104af067e45395d91f6e322",
    "norm.telomere.out": "1d6fbd74a1b405f63d342abdc0700effdc6dc9ebf892b5adcc6fea0fc83fb239",
    "norm.telomere.vcf": "8958cfe7ff4c6ebda22a0ecf88c269b5d61acff3ea1e058453081943255959a3",
    "norm.vcf": "9362223e59ed2b1621f3eb8ea332892900b0414058d8af3feb8ad672e225ddbd",
}


def fasta(name):
    sequences, current = {}, None
    for line in (FIXTURES / name).read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            sequences[current] = []
        elif current is not None:
            sequences[current].append(line.strip())
    return {key: "".join(value) for key, value in sequences.items()}


def records(name):
    rows = []
    for line in (FIXTURES / name).read_text().splitlines():
        if line and not line.startswith("#"):
            fields = line.split("\t")
            rows.append((fields[0], int(fields[1]), fields[3], fields[4]))
    return rows


def synapse_normalized(vcf, reference, *, skipped=None):
    variants, _, _ = V.parse_vcf((FIXTURES / vcf).read_text(), skipped=skipped)
    normalized = V.normalize_variants(variants, contig_references=fasta(reference))
    return [(v.chromosome, v.position, v.reference, v.alternative) for v in normalized]


class FixtureIntegrityTests(unittest.TestCase):
    def test_fixtures_are_the_published_files(self):
        for name, digest in SHA256.items():
            self.assertEqual(hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest(), digest, name)


class BcftoolsNormComparisonTests(unittest.TestCase):
    def test_split_then_normalize_matches_bcftools(self):
        """`bcftools norm -m- -f norm.fa norm.split.vcf`: 29 of 29 alleles identical.

        Covers left-shifting through homopolymer and dinucleotide repeats,
        trimming of MNPs to SNPs, an insertion at position 1 (anchored on the
        following base), a 200-nt deletion and repeated records.
        """
        skipped = []
        synapse = Counter(synapse_normalized("norm.split.vcf", "norm.fa", skipped=skipped))
        bcftools = Counter(records("norm.split.and.norm.out"))
        # bcftools keeps the reference-only record 20:59 AG>. (trimmed to A>.);
        # Synapse does not count it as a variant and reports it as skipped.
        self.assertEqual(bcftools - synapse, Counter({("20", 59, "A", "."): 1}))
        self.assertEqual(synapse - bcftools, Counter())
        self.assertEqual(sum(synapse.values()), 29)
        self.assertEqual(skipped, ["20:59 AG>."])

    def test_telomeric_repeat_shifts_to_position_one(self):
        self.assertEqual(synapse_normalized("norm.telomere.vcf", "norm.fa"), records("norm.telomere.out"))

    def test_already_left_aligned_records_are_unchanged(self):
        self.assertEqual(synapse_normalized("norm.right-align.vcf", "norm.right-align.fa"),
                         records("norm.right-align.1.out"))

    def test_ref_check_flags_exactly_the_record_bcftools_excludes(self):
        """`bcftools norm -cx` drops 1:60 A>T (REF mismatch); Synapse flags the same record."""
        variants, _, _ = V.parse_vcf((FIXTURES / "norm.vcf").read_text(), split_multiallelic=False)
        check = V.check_reference_alleles(variants, fasta("norm.fa"))
        self.assertEqual(check["mismatches"], [("1", 60, "A", "T")])
        self.assertNotIn(("1", 60), {(chrom, pos) for chrom, pos, _, _ in records("norm.out")})

    def test_api_reports_reference_only_records(self):
        from core import api
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "norm.split.vcf"
            path.write_bytes((FIXTURES / "norm.split.vcf").read_bytes())
            result = api.analyze_vcf(str(path), contig_references=fasta("norm.fa"))
        self.assertEqual(result.data["reference_only_records"], ["20:59 AG>."])
        self.assertIn("1 record(s) without an ALT allele (ALT '.', reference-only sites) are not "
                      "variants and were not counted.", result.warnings)
        self.assertEqual(len(result.data["variants"]), 29)

    def test_variant_menu_reports_reference_only_records(self):
        import builtins
        import contextlib
        import io
        from unittest import mock
        from tests.test_terminal_ui import _scripted
        from ui import variant_analysis_menu
        output = io.StringIO()
        path = str(FIXTURES / "norm.split.vcf")
        with mock.patch.object(builtins, "input", _scripted(["1", path, "0"])), \
                contextlib.redirect_stdout(output):
            variant_analysis_menu.run_variant_analysis_ui()
        self.assertIn("Loaded 29 variant(s).", output.getvalue())
        self.assertIn("NOTE: 1 record(s) without an ALT allele (ALT '.', reference-only sites) "
                      "are not variants and were not loaded.", output.getvalue())


if __name__ == "__main__":
    unittest.main()
