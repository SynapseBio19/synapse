# data/restriction_enzymes.py


# ============================================================
# RESTRICTION ENZYME DATABASE
# ============================================================
#
# Coordinate convention:
#   cut_top    = number of bases of the recognition site 5' of the
#                top-strand cut (G^AATTC -> 1).
#   cut_bottom = same position measured on the top-strand coordinate for
#                the bottom-strand cut (G^AATTC -> 5).
#   cut_top < cut_bottom  -> 5' overhang; cut_top > cut_bottom -> 3' overhang.
#
# Source and review: recognition sites and cut positions follow REBASE /
# manufacturer documentation for these commercially standard enzymes.
# Reviewed 2026-10-29 (Synapse audit): ClaI corrected from A^TCGAT to
# AT^CGAT (2-nt 5' overhang CG). All other entries unchanged.
#
# Methylation is not modelled by the digest engine. Entries with a known
# methylation dependence carry a "methylation_note" for display.

ENZYME_DATABASE_VERSION = "synapse-2026-10-29"

RESTRICTION_ENZYMES = {

    # ========================================================
    # 6-BP RECOGNITION SITES
    # ========================================================

    "EcoRI": {
        "name": "EcoRI",
        "recognition_sequence": "GAATTC",
        "cut_top": 1,
        "cut_bottom": 5,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Escherichia coli",
    },

    "BamHI": {
        "name": "BamHI",
        "recognition_sequence": "GGATCC",
        "cut_top": 1,
        "cut_bottom": 5,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Bacillus amyloliquefaciens",
    },

    "HindIII": {
        "name": "HindIII",
        "recognition_sequence": "AAGCTT",
        "cut_top": 1,
        "cut_bottom": 5,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Haemophilus influenzae",
    },

    "XhoI": {
        "name": "XhoI",
        "recognition_sequence": "CTCGAG",
        "cut_top": 1,
        "cut_bottom": 5,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Xanthomonas holcicola",
    },

    "SalI": {
        "name": "SalI",
        "recognition_sequence": "GTCGAC",
        "cut_top": 1,
        "cut_bottom": 5,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Streptomyces albus",
    },

    "PstI": {
        "name": "PstI",
        "recognition_sequence": "CTGCAG",
        "cut_top": 5,
        "cut_bottom": 1,
        "overhang_type": "3_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Providencia stuartii",
    },

    "KpnI": {
        "name": "KpnI",
        "recognition_sequence": "GGTACC",
        "cut_top": 5,
        "cut_bottom": 1,
        "overhang_type": "3_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Klebsiella pneumoniae",
    },

    "SacI": {
        "name": "SacI",
        "recognition_sequence": "GAGCTC",
        "cut_top": 5,
        "cut_bottom": 1,
        "overhang_type": "3_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Streptomyces achromogenes",
    },

    "SphI": {
        "name": "SphI",
        "recognition_sequence": "GCATGC",
        "cut_top": 5,
        "cut_bottom": 1,
        "overhang_type": "3_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Streptomyces phaeochromogenes",
    },

    "NdeI": {
        "name": "NdeI",
        "recognition_sequence": "CATATG",
        "cut_top": 2,
        "cut_bottom": 4,
        "overhang_type": "5_prime",
        "overhang_length": 2,
        "is_blunt": False,
        "source_organism": "Neisseria denitrificans",
    },

    "NcoI": {
        "name": "NcoI",
        "recognition_sequence": "CCATGG",
        "cut_top": 1,
        "cut_bottom": 5,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Neisseria corrodens",
    },

    "XbaI": {
        "name": "XbaI",
        "recognition_sequence": "TCTAGA",
        "methylation_note": "Blocked by overlapping Dam methylation (TCTAGATC / GATCTAGA).",
        "cut_top": 1,
        "cut_bottom": 5,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Xanthomonas badrii",
    },

    "SpeI": {
        "name": "SpeI",
        "recognition_sequence": "ACTAGT",
        "cut_top": 1,
        "cut_bottom": 5,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "source_organism": "Sphaerotilus natans",
    },

    "ClaI": {
        "name": "ClaI",
        "recognition_sequence": "ATCGAT",
        "cut_top": 2,
        "cut_bottom": 4,
        "overhang_type": "5_prime",
        "overhang_length": 2,
        "methylation_note": "Blocked by overlapping Dam methylation (GATCGAT / ATCGATC).",
        "is_blunt": False,
        "source_organism": "Caryophanon latum",
    },

    "EcoRV": {
        "name": "EcoRV",
        "recognition_sequence": "GATATC",
        "cut_top": 3,
        "cut_bottom": 3,
        "overhang_type": "blunt",
        "overhang_length": 0,
        "is_blunt": True,
        "source_organism": "Escherichia coli",
    },

    "SmaI": {
        "name": "SmaI",
        "recognition_sequence": "CCCGGG",
        "cut_top": 3,
        "cut_bottom": 3,
        "overhang_type": "blunt",
        "overhang_length": 0,
        "is_blunt": True,
        "source_organism": "Serratia marcescens",
    },

    "AluI": {
        "name": "AluI",
        "recognition_sequence": "AGCT",
        "cut_top": 2,
        "cut_bottom": 2,
        "overhang_type": "blunt",
        "overhang_length": 0,
        "is_blunt": True,
        "source_organism": "Arthrobacter luteus",
    },

    "HaeIII": {
        "name": "HaeIII",
        "recognition_sequence": "GGCC",
        "cut_top": 2,
        "cut_bottom": 2,
        "overhang_type": "blunt",
        "overhang_length": 0,
        "is_blunt": True,
        "source_organism": "Haemophilus aegyptius",
    },

    # ========================================================
    # FREQUENTLY USED 4-BP ENZYMES
    # ========================================================

    "MspI": {
        "name": "MspI",
        "recognition_sequence": "CCGG",
        "cut_top": 1,
        "cut_bottom": 3,
        "overhang_type": "5_prime",
        "overhang_length": 2,
        "is_blunt": False,
        "source_organism": "Moraxella species",
    },

    "HpaII": {
        "name": "HpaII",
        "recognition_sequence": "CCGG",
        "cut_top": 1,
        "cut_bottom": 3,
        "overhang_type": "5_prime",
        "overhang_length": 2,
        "is_blunt": False,
        "source_organism": "Haemophilus parainfluenzae",
        "methylation_note": "Blocked by CpG methylation (CmCGG); MspI is the insensitive isoschizomer.",
    },

    "TaqI": {
        "name": "TaqI",
        "recognition_sequence": "TCGA",
        "cut_top": 1,
        "cut_bottom": 3,
        "overhang_type": "5_prime",
        "overhang_length": 2,
        "is_blunt": False,
        "source_organism": "Thermus aquaticus",
    },

    "DpnI": {
        "name": "DpnI",
        "recognition_sequence": "GATC",
        "cut_top": 2,
        "cut_bottom": 2,
        "overhang_type": "blunt",
        "overhang_length": 0,
        "is_blunt": True,
        "source_organism": "Diplococcus pneumoniae",
        "methylation_note": "Cleaves only Dam-methylated GATC (Gm6ATC); unmethylated PCR products are not cut.",
    },
}


# ============================================================
# TYPE IIS ENZYMES (audit P2-7)
# ============================================================
#
# A Type IIS enzyme cuts outside its asymmetric recognition site. The site
# therefore occurs in two orientations, and the overhang is whatever
# sequence lies between the two cuts. cut_top / cut_bottom follow the
# convention above, measured from the first base of the site in its
# forward orientation, so they exceed the site length
# (BsaI GGTCTC(1/5): 6 + 1 and 6 + 5).
#
# The entries are kept apart from RESTRICTION_ENZYMES on purpose: the
# cloning engine gives every enzyme one fixed overhang, which a Type IIS
# enzyme does not have. They take part in digests (both strands) but not
# in the cloning compatibility model.
#
# Source: the REBASE-derived enzyme data shipped with Biopython 1.88
# (Bio.Restriction, elucidate()): BsaI GGTCTCN^NNNN_N, BsmBI
# CGTCTCN^NNNN_N, BbsI GAAGACNN^NNNN_N, SapI GCTCTTCN^NNN_N.
# tests/test_external_reference.py compares the entries and the cut
# positions on both strands with Biopython.

TYPE_IIS_DATABASE_VERSION = "synapse-2026-11-11"

TYPE_IIS_ENZYMES = {
    "BsaI": {
        "name": "BsaI",
        "recognition_sequence": "GGTCTC",
        "cut_top": 7,
        "cut_bottom": 11,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "type": "IIS",
    },

    "BsmBI": {
        "name": "BsmBI",
        "recognition_sequence": "CGTCTC",
        "cut_top": 7,
        "cut_bottom": 11,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "type": "IIS",
    },

    "BbsI": {
        "name": "BbsI",
        "recognition_sequence": "GAAGAC",
        "cut_top": 8,
        "cut_bottom": 12,
        "overhang_type": "5_prime",
        "overhang_length": 4,
        "is_blunt": False,
        "type": "IIS",
    },

    "SapI": {
        "name": "SapI",
        "recognition_sequence": "GCTCTTC",
        "cut_top": 8,
        "cut_bottom": 11,
        "overhang_type": "5_prime",
        "overhang_length": 3,
        "is_blunt": False,
        "type": "IIS",
    },
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_restriction_enzyme(
    enzyme_name
):

    return RESTRICTION_ENZYMES.get(
        enzyme_name
    )


def get_type_iis_enzyme(
    enzyme_name
):

    return TYPE_IIS_ENZYMES.get(
        enzyme_name
    )


def get_all_enzyme_names():

    return sorted(
        RESTRICTION_ENZYMES.keys()
    )


def get_blunt_end_enzymes():

    return [
        name
        for name, enzyme
        in RESTRICTION_ENZYMES.items()
        if enzyme["is_blunt"]
    ]


def get_sticky_end_enzymes():

    return [
        name
        for name, enzyme
        in RESTRICTION_ENZYMES.items()
        if not enzyme["is_blunt"]
    ]


def get_five_prime_overhang_enzymes():

    return [
        name
        for name, enzyme
        in RESTRICTION_ENZYMES.items()
        if enzyme["overhang_type"]
        == "5_prime"
    ]


def get_three_prime_overhang_enzymes():

    return [
        name
        for name, enzyme
        in RESTRICTION_ENZYMES.items()
        if enzyme["overhang_type"]
        == "3_prime"
    ]