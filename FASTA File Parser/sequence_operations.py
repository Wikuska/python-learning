from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re

def find_sequence(filename, seq_id):
    for seq_record in SeqIO.parse(filename, "fasta"):
        if seq_record.id == seq_id:
            return seq_record

def seq_length(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    length = len(seq_record)
    return length

def nucleotides_count(filename, seq_id, seq_type):
    count_dict = {}
    seq_record = find_sequence(filename, seq_id)
    if seq_type == "DNA":
        count_dict = {
            "A": seq_record.count("A"),
            "T": seq_record.count("T"),
            "C": seq_record.count("C"),
            "G": seq_record.count("G")
        }
    elif seq_type == "RNA":
        count_dict = {
            "A": seq_record.count("A"),
            "U": seq_record.count("U"),
            "C": seq_record.count("C"),
            "G": seq_record.count("G")
        }
    return count_dict
        
def gc_content(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    content = round(100 * (gc_fraction(seq_record)))
    return content

def get_complementary(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    complementary_seq = (seq_record.seq).complement()
    return str(complementary_seq)

def get_reverse_complementary(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    reverse_complementary_seq = (seq_record.seq).reverse_complement()
    return str(reverse_complementary_seq)

def transcribe_seq(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    m_rna = (seq_record.seq).transcribe()
    return str(m_rna)