from Bio.Seq import Seq

def app():
    sequence = input("Type in your DNA sequence: ").strip().upper()
    if not all(base in "ATCG" for base in sequence):
        print("Invalid sequence\n")
        app()

    seq = Seq(sequence)
    while True:
        print("\nHere are operations you can make:")
        print("1. Get complementary sequence")
        operation = int(input("What would you like to do? Choose number from options: "))
        print("--------------------------------------------------------")

        if operation == 1:
            print(seq.complement())

        print("--------------------------------------------------------")    
        go_on = input("\nWould you like to make more operations?(y/n): ").upper()
        if go_on == "N":
            break

app()
