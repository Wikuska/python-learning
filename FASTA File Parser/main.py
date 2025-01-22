from Bio import SeqIO
from operations import records_num, records_lenght, find_longest_shortest

OPERATIONS_LIST = ["Show operations list", "Get records number", "Get each record length", "Get the longest and the shortest record"]

def show_operations_list():
    print("Operations list:")
    for operation in OPERATIONS_LIST:
        print(str(OPERATIONS_LIST.index(operation) + 1) + ". " + operation)

def choose_operation():
    while True:
        try:
            operation = int(input("What would you like to do with your data? Enter operation number: "))
            if operation > len(OPERATIONS_LIST) or operation < 1:
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

    show_operations_list()
    
    while True:
        operation = choose_operation()

        if operation == 1:
            show_operations_list()

        if operation == 2:
            records = records_num(filename)
            print(f"There are {records} records in this file")

        if operation == 3:
            length_dict = records_lenght(filename)
            for key, value in length_dict.items():
                print(f"\nID: {value['id']}\nDescription: {value['description']}\nSequence length: {value["length"]}")

        if operation == 4:
            longest_shortest_dict = find_longest_shortest(filename)
            if longest_shortest_dict:
                print(f"\nLongest:\nID{longest_shortest_dict["longest"]["id"]}\nDescription: {longest_shortest_dict["longest"]["description"]}\nLength: {longest_shortest_dict["longest"]["length"]}")
                print(f"\nLongest:\nID{longest_shortest_dict["shortest"]["id"]}\nDescription: {longest_shortest_dict["shortest"]["description"]}\nLength: {longest_shortest_dict["shortest"]["length"]}")
            else:
                print("\nThis file contains only one sequence")

        while True:
            go_on = input("\nWould you like to do something else? (y/n): ").lower().strip()
            if go_on == "n" or go_on == "y":
                break
            else:
                print("Invalid input")

        if go_on == "n":
            break
        

app()
