from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import gc_fraction

def find_sequence(filename, seq_id):
    for seq_record in SeqIO.parse(filename, "fasta"):
        if seq_record.id == seq_id:
            return seq_record
        
def seq_id(id):
    return f"Sequence file id: {id}"

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
    return f"Complementary sequence:\n{str(complementary_seq)}"

def get_reverse_complementary(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    reverse_complementary_seq = (seq_record.seq).reverse_complement()
    return f"Reverse complementary sequence:\n{str(reverse_complementary_seq)}"

def transcribe_to_mrna(filename, seq_id):
    seq_record = find_sequence(filename, seq_id)
    m_rna = (seq_record.seq).transcribe()
    return f"Messenger RNA sequence:\n{str(m_rna)}"

# -----------------------------------PROTEIN ------------------------------------- #

def get_molecular_weight(filename, seq_id):
    sequence = ProteinAnalysis(str((find_sequence(filename, seq_id).seq)))
    sequence.molecular_weight()
    return f"Protein molecular weight: {sequence.molecular_weight()/1000:.2f} kDa"

def get_amino_acids_occurrence(filename, seq_id):
    sequence = ProteinAnalysis(str((find_sequence(filename, seq_id).seq)))
    acids_dict = sequence.count_amino_acids()
    perc_dict = sequence.amino_acids_percent
    output = "Amino acids occurrance:"
    for (name, occurrence), occ_perc in zip(acids_dict.items(), perc_dict.values()):
        output += f"\n{name} - {occurrence} - {occ_perc:.2f}%"
    return output

def find_isoelectric_point(filename, seq_id):
    sequence = ProteinAnalysis(str((find_sequence(filename, seq_id).seq)))
    return f"Isoelectric point: {sequence.isoelectric_point():.2f} pH"

def check_stability(filename, seq_id):
    sequence = ProteinAnalysis(str((find_sequence(filename, seq_id).seq)))
    if sequence.instability_index() > 40:
        return "Protein is unstable"
    else:
        return "Protein is stable"


