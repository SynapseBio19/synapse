dna = input("Enter DNA Sequence: ")

sequence_is_valid = True

for nucleotide in dna:
    if nucleotide not in "ATCG":
        sequence_is_valid = False
        break

if sequence_is_valid:
    a_count = dna.count("A")
    t_count = dna.count("T")
    c_count = dna.count("C")
    g_count = dna.count("G")
    sequence_length = len(dna)

    gc_content = (c_count + g_count) * 100 / sequence_length

    print("DNA Sequence:", dna)
    print("Length:", sequence_length)

    print("A count:", a_count)
    print("T count:", t_count)
    print("C count:", c_count)
    print("G count:", g_count)

    print(f"GC Content: {gc_content:.2f}%")
    print("The DNA sequence is valid.")
    
else:
    print("The DNA sequence is invalid.")

input("Press Enter to exit...")