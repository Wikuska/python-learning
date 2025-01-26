from Bio import SeqIO

def records_num(filename):
    records = 0
    for seq_record in SeqIO.parse(filename, "fasta"):
        records += 1
    return f"Records found in file: {records}"

def records_length(filename):
    records_length_output = ""
    for seq_record in SeqIO.parse(filename, "fasta"):
        records_length_output += f"{seq_record.id} - {len(seq_record)}\n"
    return records_length_output

def find_longest_shortest(filename):
    length_list = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        length_list.append(len(seq_record))
    
    if len(length_list) == 1:
        return "There is only one record in file"
    else:
        longest_id = []
        shortest_id = []

        for seq_record in SeqIO.parse(filename, "fasta"):
            if len(seq_record) == max(length_list):
                longest_id.append(seq_record.id)
            elif len(seq_record) == min(length_list):
                shortest_id.append(seq_record.id)

        longest_output = "Longest:"
        shortest_output = "Shortest:"

        for id in longest_id:
            longest_output += f"\n{id} -- {max(length_list)}"

        for id in shortest_id:
            shortest_output += f"\n{id} -- {min(length_list)}"

        return f"{longest_output}\n{shortest_output}"
        
    
    
