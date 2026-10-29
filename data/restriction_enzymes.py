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
# HELPER FUNCTIONS
# ============================================================

def get_restriction_enzyme(
    enzyme_name
):

    return RESTRICTION_ENZYMES.get(
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