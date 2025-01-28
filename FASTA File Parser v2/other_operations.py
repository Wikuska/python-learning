from Bio import SeqIO
from PySide6.QtWidgets import  QApplication, QMainWindow, QFileDialog

def check_file(app):
    file_path, _ = QFileDialog.getOpenFileName(app, "Choose a file", "", "FASTA Files (*.fasta *.fa);;All Files (*)")
    if file_path:
        try:
            SeqIO.parse(file_path, "fasta")
        except Exception:
            return None
        else:
            return file_path
        
def records_ids(filename):
    ids_list = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        ids_list.append(seq_record.id)
    return ids_list
        
def check_seq_type(filename, seq_id):
    for seq_record in SeqIO.parse(filename, "fasta"):
        if seq_record.id == seq_id:
            sequence = str(seq_record.seq).upper()
            if set(sequence).issubset(set("ATCGN")):
                return "DNA"
            elif set(sequence).issubset(set("AUCGN")):
                return "RNA"
            elif all(char.isalpha() for char in sequence):
                return "Protein"
            else:
                return "Unknown"

def seperate_longer_str(text):
    formatted_string = " ".join(text[i:i+65] for i in range(0, len(text), 65))
    return formatted_string