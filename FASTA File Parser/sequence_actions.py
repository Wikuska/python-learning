from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from other_operations import seperate_longer_str

def find_sequence(filename, seq_id):
    for seq_record in SeqIO.parse(filename, "fasta"):
        if seq_record.id == seq_id:
            return seq_record

def seq_length(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    length = len(seq_record)
    return f"Sequence length: {length}"

def nucleotides_count(filename, seq_id, seq_type):
    seq_record = find_sequence(filename, seq_id)
    a_count = seq_record.count("A")
    t_count = seq_record.count("T")
    c_count = seq_record.count("C")
    g_count = seq_record.count("G")
    u_count = seq_record.count("U")
    if seq_type == "DNA":
        return f"Nucleotides in this sequence:\nA = {a_count}\nT = {t_count}\nC = {c_count}\nG = {g_count}"
    elif seq_type == "RNA":
        return f"Nucleotides in this sequence:\nA = {a_count}\nU = {u_count}\nC = {c_count}\nG = {g_count}"
    
def gc_content(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    content = round(100 * (gc_fraction(seq_record)), 2)
    return f"GC content in this sequence: {content} %"

def get_complementary(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    complementary_seq = (seq_record.seq).complement()
    formatted_sequence = seperate_longer_str(str(complementary_seq))
    return f"Complementary sequence:\n{formatted_sequence}"

def get_reverse_complementary(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    reverse_complementary_seq = (seq_record.seq).reverse_complement()
    formatted_sequence = seperate_longer_str(str(reverse_complementary_seq))
    return f"Reverse complementary sequence:\n{formatted_sequence}"

def transcribe_to_mrna(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    m_rna = (seq_record.seq).transcribe()
    formatted_sequence = seperate_longer_str(str(m_rna))
    return f"Messenger RNA sequence:\n{formatted_sequence}"