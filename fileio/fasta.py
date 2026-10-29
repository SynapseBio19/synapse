import os


def load_fasta(filename):

    sequences = []
    header = None
    dna_parts = []

    with open(filename, "r") as file:

        for line_number, line in enumerate(file, start=1):

            line = line.strip()

            # Ignore empty lines
            if not line:
                continue

            # FASTA header
            if line.startswith(">"):

                new_header = line[1:].strip()

                # Header must not be empty
                if not new_header:
                    raise ValueError(
                        f"FASTA format error at line {line_number}: "
                        "empty header."
                    )

                # Save previous record
                if header is not None:

                    if not dna_parts:
                        raise ValueError(
                            f"FASTA format error: "
                            f"'{header}' contains no sequence."
                        )

                    sequences.append({
                        "header": header,
                        "sequence": "".join(dna_parts)
                    })

                header = new_header
                dna_parts = []

            else:

                # Sequence before first header
                if header is None:
                    raise ValueError(
                        f"FASTA format error at line {line_number}: "
                        "sequence found before header."
                    )

                dna_parts.append(line)

    # No header was ever found
    if header is None:
        raise ValueError(
            "FASTA file contains no sequences."
        )

    # Last header has no sequence
    if not dna_parts:
        raise ValueError(
            f"FASTA format error: "
            f"'{header}' contains no sequence."
        )

    # Save final record
    sequences.append({
        "header": header,
        "sequence": "".join(dna_parts)
    })

    return sequences