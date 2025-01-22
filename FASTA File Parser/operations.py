from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re

def records_num(filename):
    records = 0
    for seq_record in SeqIO.parse(filename, "fasta"):
        records += 1
    return records

def records_lenght(filename):
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
    records_length_dict = records_lenght(filename)
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

