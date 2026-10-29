def show_menu():

    print()
    print("=" * 60)
    print("Choose Analysis")
    print("=" * 60)

    print("[1] Sequence Overview")
    print("[2] ORF Analysis")
    print("[3] RNA Analysis")
    print("[4] Protein Analysis")

    print("-" * 60)

    print("[5] Restriction & Cloning Analysis")
    print("[6] Primer Design")
    print("[7] Mutation & Variant Analysis")
    print("[8] CRISPR Guide Design")

    print("-" * 60)

    print("[9] BLAST Search")
    print("[10] Literature Search")
    print("[11] AI Experiment Planner")

    print("-" * 60)

    print("[S] Select Another Sequence")

    print("-" * 60)

    print("[0] Exit")

    return input("\nSelection: ")


def select_gene(sequences):

    print("=" * 60)
    print("SYNAPSE")
    print("=" * 60)

    print("\nLoaded Genes\n")

    for index, record in enumerate(sequences, start=1):
        print(f"[{index}] {record['header']}")

    while True:

        choice = input("\nSelect a gene: ")

        if choice.isdigit():

            choice = int(choice)

            if 1 <= choice <= len(sequences):
                return choice, sequences[choice - 1]

        print("Invalid selection.")