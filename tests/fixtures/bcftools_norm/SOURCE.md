# bcftools norm test vectors

Verbatim copies of input and expected-output files from the bcftools test
suite, used by `tests/test_bcftools_norm_reference.py` to compare Synapse's
variant normalization with the documented output of `bcftools norm`.

- Source: <https://github.com/samtools/bcftools/tree/1.25/test>
  (release tag `1.25`, published 2026-10-27), raw files from
  `https://raw.githubusercontent.com/samtools/bcftools/1.25/test/<name>`.
- Downloaded: 2026-11-13.
- Licence: bcftools is available under the MIT/Expat licence or the GPL
  (Copyright (C) 2012-2024 Genome Research Ltd.); these test files are
  redistributed here under the MIT/Expat terms.
- How bcftools produced each output (from `test/test.pl` of the same tag):

| Expected output | Input | Reference | `bcftools norm` arguments |
|---|---|---|---|
| `norm.out` | `norm.vcf` | `norm.fa` | `-cx` (check REF, exclude mismatches; multi-allelic records normalized jointly) |
| `norm.split.and.norm.out` | `norm.split.vcf` | `norm.fa` | `-m-` (split multi-allelic records, then normalize) |
| `norm.telomere.out` | `norm.telomere.vcf` | `norm.fa` | none (default left-alignment) |
| `norm.right-align.1.out` | `norm.right-align.vcf` | `norm.right-align.fa` | `--old-rec-tag ORI` (default left-alignment) |

SHA-256 (checked by the test):

```
da20982fbeb50362b580d8d2e06c492bfe221cb92020729a3e3071147764b060  norm.fa
dba0b08005378b011f90695ee1ac5786eaf9b0caab767c6f7e28a54fd87a6b0a  norm.out
f16444cd7521694a32e8e18ce041afc63301800e200f7892d8def419c7722bc6  norm.right-align.1.out
0eb5e4aefeabacac15326d354dbf4a478415a92f31d926332457e8048547d497  norm.right-align.fa
85004daeb31afd43dba635150c8cea707e4963a21ae5f194a90e3018ff015a59  norm.right-align.vcf
0e4412531059a5e49757a27a9d96e00dc33375301dcb4eb71b1f34e11110b221  norm.split.and.norm.out
0363806f590a1ce92c7e2a8c3e16d976577e91fec104af067e45395d91f6e322  norm.split.vcf
1d6fbd74a1b405f63d342abdc0700effdc6dc9ebf892b5adcc6fea0fc83fb239  norm.telomere.out
8958cfe7ff4c6ebda22a0ecf88c269b5d61acff3ea1e058453081943255959a3  norm.telomere.vcf
9362223e59ed2b1621f3eb8ea332892900b0414058d8af3feb8ad672e225ddbd  norm.vcf
```

bcftools itself was not run (no native Windows build; WSL is not
installed). The comparison is against the outputs that the bcftools
project publishes and checks in its own test suite.
