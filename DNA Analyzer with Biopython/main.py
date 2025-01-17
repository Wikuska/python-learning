from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

def app():
    sequence = input("Type in your DNA sequence: ").strip().upper()
    if not all(base in "ATCG" for base in sequence):
        print("Invalid sequence\n")
        app()

    seq = Seq(sequence)
    while True:
        print("\nHere are operations you can make:")
        print("1. Get complementary sequence")
        print("2. Get sequence length")
        print("3. Get amount of each nucleotides")
        print("4. Get GC %")
        operation = int(input("What would you like to do? Choose number from options: "))
        print("--------------------------------------------------------")

        if operation == 1:
            print(seq.complement())
        elif operation == 2:
            print(len(seq))
        elif operation == 3:
            print(f"A: {seq.count("A")}\nT: {seq.count("T")}\nC: {seq.count("C")}\nG: {seq.count("G")}")
        elif operation == 4:
            print(round(100 * (gc_fraction(seq)), 1))

        print("--------------------------------------------------------")    
        go_on = input("\nWould you like to make more operations?(y/n): ").upper()
        if go_on == "N":
            break

app()
