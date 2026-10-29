"""Standard genetic code (NCBI table 1), RNA alphabet; derived from core.seqtools."""

from core.seqtools import codon_dict

CODON_TABLE = codon_dict(1, rna=True, stop_symbol="STOP")