from ui import style


def show_menu():

    style.header("Main menu")
    style.menu([
        ("Sequence analysis", [
            (1, "Sequence Overview"),
            (2, "ORF Analysis"),
            (3, "RNA Analysis"),
            (4, "Protein Analysis"),
        ]),
        ("Design and comparison", [
            (5, "Restriction & Cloning Analysis"),
            (6, "Primer Design"),
            (7, "Mutation & Variant Analysis"),
            (8, "CRISPR Guide Design"),
        ]),
        ("Search and planning", [
            (9, "BLAST Search"),
            (10, "Literature Search"),
            (11, "AI Experiment Planner (not implemented)"),
        ]),
        ("Session", [
            ("P", "Pipeline Status & Stored Results"),
            ("R", "Export Analysis Report (JSON / Markdown)"),
            ("S", "Select Another Sequence"),
        ]),
    ], zero=style.EXIT_LABEL)

    return style.ask_choice()


def select_gene(sequences):

    style.header("Loaded sequences")

    entries = []
    for index, record in enumerate(sequences, start=1):
        topology = record.get("topology", "linear")
        length = len(record.get("sequence", ""))
        unit = "aa" if record.get("alphabet") == "PROTEIN" else "bp"
        entries.append((index, f"{record['header']}  ({length} {unit}, {topology})"))

    # The selector cannot be left with 0: a sequence or [F] must be chosen.
    style.menu([("Sequences", entries), ("File", [("F", "Load another file")])], zero=None)

    while True:

        choice = style.ask_choice()

        if choice.upper() == "F":
            return None, None

        if choice.isdigit():

            choice = int(choice)

            if 1 <= choice <= len(sequences):
                return choice, sequences[choice - 1]

        style.invalid_selection()
