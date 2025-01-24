from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re

def records_num(filename):
    records = 0
    for seq_record in SeqIO.parse(filename, "fasta"):
        records += 1
    return records

def records_length(filename):
    recods_length_dict = {}
    i = 1
    for seq_record in SeqIO.parse(filename, "fasta"):
        sequence_id = seq_record.id
        description = seq_record.description
        length = len(seq_record)
        recods_length_dict[i] = {
            "id": sequence_id,
            "description" : description,
            "length": length
        }
        i += 1
    return recods_length_dict

def find_longest_shortest(filename):
    records_length_dict = records_length(filename)
    if len(records_length_dict) < 2:
        return False
    else:
        longest_shortest_dict = {}
        for key, value in records_length_dict.items():
            if key == 1:
                longest_shortest_dict["longest"] = {"id": value["id"], "description": value["description"], "length": value["length"]}
                longest_shortest_dict["shortest"] = {"id": value["id"], "description": value["description"], "length": value["length"]}
            else:
                if value["length"] > longest_shortest_dict["longest"]["length"]:
                    longest_shortest_dict["longest"]["id"] = value["id"]
                    longest_shortest_dict["longest"]["description"] = value["description"]
                    longest_shortest_dict["longest"]["length"] = value["length"]
                elif value["length"] < longest_shortest_dict["shortest"]["length"]:
                    longest_shortest_dict["shortest"]["id"] = value["id"]
                    longest_shortest_dict["shortest"]["description"] = value["description"]
                    longest_shortest_dict["shortest"]["length"] = value["length"]
        
        return longest_shortest_dict
    
def choose_record(filename):
    print("Here is list of records from this file you can work on:")
    record_list = []
    i = 1
    for seq_record in SeqIO.parse(filename, "fasta"):
        print(f"{i}: {seq_record.description}")
        record_list.append(str(seq_record.seq))
        i += 1

    while True:
        try:
            sequence_num = int(input("Enter number of sequence you want to use: ").strip())
            if sequence_num > len(record_list) or sequence_num < 1:
                raise IndexError
        except ValueError:
            print("Please enter a number.")
        except IndexError:
            print("No such operation avaliable.")
        except Exception as e:
            print(f"An unexpected error occured: {e}")
        else:
            sequence_string = record_list[sequence_num - 1]
            return sequence_string
