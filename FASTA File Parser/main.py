from Bio import SeqIO

OPERATIONS_LIST = ["Get records number", "Get each record length"]

def app():
    print("Welcome in FASTA File Parser")

    while True:
        filename = input("Please enter name of the file you want to work on (including file extension): ").strip()
        try:
            SeqIO.parse(filename, "fasta")
            print(f"File '{filename}' successfully found and loaded.")
            break
        except FileNotFoundError:
            print("File like this was not found. Make sure you provided correct name including extension.")
        except Exception as e:
            print(f"An unexpected error occured: {e}")

    print("Operations list:")
    for operation in OPERATIONS_LIST:
        print(str(OPERATIONS_LIST.index(operation) + 1) + ". " + operation)

app()
