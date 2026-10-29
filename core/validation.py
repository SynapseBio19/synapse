def normalize_dna(dna):

    dna = dna.upper()
    dna = dna.replace(" ", "")
    dna = dna.replace("\n", "")
    dna = dna.replace("\t", "")

    return dna

def validate_dna(dna):

    valid_nucleotides = "ATCGRYSWKMBDHVN"

    for nucleotide in dna:

        if nucleotide not in valid_nucleotides:
            return False

    return True

def get_dna_errors(dna):

    valid_nucleotides = "ATCGRYSWKMBDHVN"

    errors = []

    for position, nucleotide in enumerate(dna, start=1):

        if nucleotide not in valid_nucleotides:

            errors.append({
                "character": nucleotide,
                "position": position
            })

    return errors