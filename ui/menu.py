from ui.style import major_rule, minor_rule
def show_menu():

    print()
    print(major_rule())
    print("Choose Analysis")
    print(major_rule())

    print("[1] Sequence Overview")
    print("[2] ORF Analysis")
    print("[3] RNA Analysis")
    print("[4] Protein Analysis")

    print(minor_rule())

    print("[5] Restriction & Cloning Analysis")
    print("[6] Primer Design")
    print("[7] Mutation & Variant Analysis")
    print("[8] CRISPR Guide Design")

    print(minor_rule())

    print("[9] BLAST Search")
    print("[10] Literature Search")
    print("[11] AI Experiment Planner (not implemented)")

    print(minor_rule())

    print("[P] Pipeline Status & Stored Results")
    print("[R] Export Analysis Report (JSON / Markdown)")
    print("[S] Select Another Sequence")

    print(minor_rule())

    print("[0] Exit")

    return input("\nSelection: ")


def select_gene(sequences):

    print(major_rule())
    print("SYNAPSE")
    print(major_rule())

    print("\nLoaded Genes\n")

    for index, record in enumerate(sequences, start=1):
        topology = record.get("topology", "linear")
        length = len(record.get("sequence", ""))
        unit = "aa" if record.get("alphabet") == "PROTEIN" else "bp"
        print(f"[{index}] {record['header']}  ({length} {unit}, {topology})")

    print("\n[F] Load another file")

    while True:

        choice = input("\nSelect a sequence: ").strip()

        if choice.upper() == "F":
            return None, None

        if choice.isdigit():

            choice = int(choice)

            if 1 <= choice <= len(sequences):
                return choice, sequences[choice - 1]

        print("Invalid selection.")