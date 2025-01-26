from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re

        
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