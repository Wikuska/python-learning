from Bio import SeqIO
from operations import records_num, records_lenght, find_longest_shortest, choose_record

FILE_OPERATIONS_LIST = ["Show operations list", "Get records number", "Get each record length", "Get the longest and the shortest record", "Work on particular record"]
SEQUENCE_OPERATIONS_LIST = ["Show operations list", "Choose another sequence", "Go back to file operations"]

def show_operations_list(operation_class):
    if operation_class == "file":
        print("File operations list:")
        for operation in FILE_OPERATIONS_LIST:
            print(str(FILE_OPERATIONS_LIST.index(operation) + 1) + ". " + operation)
    elif operation_class == "sequence":
        print("Sequence operations list:")
        for operation in SEQUENCE_OPERATIONS_LIST:
            print(str(SEQUENCE_OPERATIONS_LIST.index(operation) + 1) + ". " + operation)

def choose_operation(operation_class):
    while True:
        try:
            if operation_class == "file":
                operation = int(input("What would you like to do with your data? Enter operation number: "))
                if operation > len(FILE_OPERATIONS_LIST) or operation < 1:
                    raise IndexError
            elif operation_class == "sequence":
                operation = int(input("What would you like to do with this sequence? Enter operation number: "))
                if operation > len(SEQUENCE_OPERATIONS_LIST) or operation < 1:
                    raise IndexError
        except ValueError:
            print("Please enter a number.")
        except IndexError:
            print("No such operation avaliable.")
        except Exception as e:
            print(f"An unexpected error occured: {e}")
        else:
            return operation


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

    show_operations_list("file")
    
    while True:
        file_operation = choose_operation("file")

        if file_operation == 1:
            show_operations_list("file")

        if file_operation == 2:
            records = records_num(filename)
            print(f"There are {records} records in this file")

        if file_operation == 3:
            length_dict = records_lenght(filename)
            for key, value in length_dict.items():
                print(f"\nID: {value['id']}\nDescription: {value['description']}\nSequence length: {value["length"]}")

        if file_operation == 4:
            longest_shortest_dict = find_longest_shortest(filename)
            if longest_shortest_dict:
                print(f"\nLongest:\nID{longest_shortest_dict["longest"]["id"]}\nDescription: {longest_shortest_dict["longest"]["description"]}\nLength: {longest_shortest_dict["longest"]["length"]}")
                print(f"\nShortest:\nID{longest_shortest_dict["shortest"]["id"]}\nDescription: {longest_shortest_dict["shortest"]["description"]}\nLength: {longest_shortest_dict["shortest"]["length"]}")
            else:
                print("\nThis file contains only one sequence")

        if file_operation == 5:
            while True:
                sequence_string = choose_record(filename)
                print(f"Chosen sequence: {sequence_string}")
                show_operations_list("sequence")
                while True:
                    sequence_operation = choose_operation("sequence")

                    if sequence_operation == 1:
                        show_operations_list("sequence")

                    if sequence_operation == 2:
                        break

                    if sequence_operation == 3:
                        break

                if sequence_operation == 3:
                    break

app()
