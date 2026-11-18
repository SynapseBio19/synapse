"""Fixed configuration of the BLAST module: NCBI endpoint and rate limits,
programs, remote databases and the local database catalog."""

from __future__ import annotations

import os
from core.alphabets import PROTEIN_IUPAC
from pathlib import Path

from analyses._blast.paths import ANALYSES_DIR
from core.paths import HOME_ENV, is_source_checkout, user_data_dir


# ============================================================================
# CONSTANTS
# ============================================================================

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
DEFAULT_TOOL = os.environ.get("SYNAPSE_NCBI_TOOL", "Synapse").strip() or "Synapse"
DEFAULT_EMAIL = os.environ.get("SYNAPSE_NCBI_EMAIL", "").strip()
API_TIMEOUT_SECONDS = 30
# NCBI asks API clients not to contact the BLAST server more often than once
# every 10 seconds and not to poll an individual RID more often than once a
# minute. These limits are enforced by the client instead of relying on callers.
NCBI_MIN_REQUEST_INTERVAL_SECONDS = 10.0
NCBI_MIN_POLL_INTERVAL_SECONDS = 60.0
DEFAULT_POLL_INTERVAL_SECONDS = NCBI_MIN_POLL_INTERVAL_SECONDS
DEFAULT_MAX_WAIT_SECONDS = 300.0
DEFAULT_MAX_HITS = 50

# Evidence bands are assigned only to hits whose best E-value is at or below
# this threshold (audit P0-1). 1e-3 is a common, conservative convention for
# inferring homology from a single search; it is NOT universally valid. The
# E-value depends on database size and query length, so the threshold is a
# parameter of calculate_evidence_score()/build_evidence_profile().
DEFAULT_EVIDENCE_EVALUE_THRESHOLD = 1e-3


# Local BLAST database storage. The path can be moved without changing Synapse
# code by setting SYNAPSE_BLAST_DB_ROOT in the environment.
LOCAL_DB_ROOT_ENV = "SYNAPSE_BLAST_DB_ROOT"
LOCAL_DB_MANIFEST_NAME = "manifest.json"


def _resolve_default_local_db_root() -> Path:
    """Resolve Synapse's BLAST database root independently of the CWD.

    Synapse may import ``blast.py`` from a subdirectory such as ``analyses``.
    Using ``Path(__file__).parent / data / blast_db`` in that situation creates
    the wrong location (for example ``analyses/data/blast_db``). We therefore
    walk upward and prefer an existing project-level ``data/blast_db``.

    The environment variable remains the highest-priority override.
    """
    configured = os.environ.get(LOCAL_DB_ROOT_ENV, "").strip()
    if configured:
        return Path(configured).expanduser().resolve()

    # An installed Synapse lives in site-packages, which also holds a "data"
    # package; databases must not land there (Phase 9). Use the per-user
    # folder instead. From the source folder nothing changes.
    if not is_source_checkout() or os.environ.get(HOME_ENV, "").strip():
        return (user_data_dir() / "blast_db").resolve()

    module_dir = ANALYSES_DIR
    candidates = [module_dir, *module_dir.parents]

    # Prefer an already existing BLAST root. This preserves databases that
    # were created before Synapse starts and avoids creating a second root.
    for candidate in candidates:
        blast_root = candidate / "data" / "blast_db"
        if blast_root.is_dir():
            return blast_root.resolve()

    # If no BLAST directory exists yet, identify the project root by common
    # Synapse markers before falling back to the module directory.
    for candidate in candidates:
        if (candidate / "main.py").is_file() or (candidate / "data").is_dir():
            return (candidate / "data" / "blast_db").resolve()

    return (module_dir / "data" / "blast_db").resolve()


DEFAULT_LOCAL_DB_ROOT = _resolve_default_local_db_root()

# NCBI preformatted database catalog. These names are passed directly to
# update_blastdb.pl. The catalog is deliberately conservative: Synapse never
# downloads a large database automatically.
LOCAL_DATABASE_CATALOG = {
    "refseq_protein": {
        "type": "protein",
        "description": "NCBI RefSeq Protein",
        "source": "NCBI preformatted BLAST database",
        "recommended": True,
    },
    "refseq_rna": {
        "type": "nucleotide",
        "description": "NCBI RefSeq RNA",
        "source": "NCBI preformatted BLAST database",
        "recommended": True,
    },
    "refseq_genomic": {
        "type": "nucleotide",
        "description": "NCBI RefSeq genomic",
        "source": "NCBI preformatted BLAST database",
        "recommended": False,
    },
    "swissprot": {
        "type": "protein",
        "description": "UniProtKB/Swiss-Prot via NCBI",
        "source": "NCBI preformatted BLAST database",
        "recommended": True,
    },
    "core_nt": {
        "type": "nucleotide",
        "description": "NCBI core nucleotide collection (NCBI's default since 2024)",
        "source": "NCBI preformatted BLAST database",
        "recommended": False,
    },
    "nt": {
        "type": "nucleotide",
        "description": "NCBI nucleotide collection",
        "source": "NCBI preformatted BLAST database",
        "recommended": False,
    },
    "nr": {
        "type": "protein",
        "description": "NCBI non-redundant protein collection",
        "source": "NCBI preformatted BLAST database",
        "recommended": False,
    },
}

BLAST_EXECUTABLES = (
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "blastdbcmd",
    "makeblastdb",
)
BLAST_UPDATE_SCRIPT = "update_blastdb.pl"

DNA_BASES = set("ACGT")
# IUPAC protein letters incl. J, U (Sec) and O (Pyl); BLAST+ 2.17.0 accepts
# all three in blastp queries (O is scored as X, with a warning).
PROTEIN_AMINO_ACIDS = set(PROTEIN_IUPAC)
AMBIGUOUS_DNA = set("NRYKMSWBDHV")

PROGRAMS = {
    "blastn": {
        "query_type": "DNA",
        "database_type": "nucleotide",
        "description": "Nucleotide query against a nucleotide database.",
    },
    "blastp": {
        "query_type": "PROTEIN",
        "database_type": "protein",
        "description": "Protein query against a protein database.",
    },
    "blastx": {
        "query_type": "DNA",
        "database_type": "protein",
        "description": "Six-frame translated nucleotide query against protein.",
    },
    "tblastn": {
        "query_type": "PROTEIN",
        "database_type": "nucleotide",
        "description": "Protein query against six-frame translated nucleotide sequences.",
    },
    "tblastx": {
        "query_type": "DNA",
        "database_type": "nucleotide",
        "description": "Translated nucleotide query against translated nucleotide sequences.",
    },
}

# Query types use the vocabulary of BlastQuery.sequence_type.
QUERY_TYPES = ("DNA", "PROTEIN")


DATABASES = {
    "core_nt": {
        "type": "nucleotide",
        "description": "NCBI core nucleotide collection (NCBI's default since 2024).",
    },
    "nt": {
        "type": "nucleotide",
        "description": "NCBI nucleotide collection.",
    },
    "refseq_rna": {
        "type": "nucleotide",
        "description": "NCBI RefSeq RNA collection.",
    },
    "refseq_genomic": {
        "type": "nucleotide",
        "description": "NCBI RefSeq genomic collection.",
    },
    "nr": {
        "type": "protein",
        "description": "NCBI non-redundant protein collection.",
    },
    "refseq_protein": {
        "type": "protein",
        "description": "NCBI RefSeq protein collection.",
    },
    "swissprot": {
        "type": "protein",
        "description": "UniProtKB/Swiss-Prot protein collection as exposed by NCBI.",
    },
}

# NCBI BLAST XML 2/legacy XML namespaces can vary. We deliberately parse
# local tag names so the parser is resilient to namespace prefixes.
BLAST_TAG_ALIASES = {
    "Iteration",
    "Hit",
    "Hsp",
    "Hit_id",
    "Hit_def",
    "Hit_accession",
    "Hit_len",
    "Hsp_num",
    "Hsp_bit-score",
    "Hsp_score",
    "Hsp_evalue",
    "Hsp_query-from",
    "Hsp_query-to",
    "Hsp_hit-from",
    "Hsp_hit-to",
    "Hsp_query-frame",
    "Hsp_hit-frame",
    "Hsp_identity",
    "Hsp_positive",
    "Hsp_gaps",
    "Hsp_align-len",
    "Hsp_qseq",
    "Hsp_hseq",
    "Hsp_midline",
    "BlastOutput_program",
    "BlastOutput_version",
    "BlastOutput_db",
    "Parameters_expect",
    "Parameters_matrix",
    "Parameters_gap-open",
    "Parameters_gap-extend",
    "Statistics_db-num",
    "Statistics_db-len",
    "Statistics_hsp-len",
}
