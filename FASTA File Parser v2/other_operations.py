from Bio import SeqIO
from PySide6.QtWidgets import QFileDialog, QMessageBox

def check_file(app, filename_label):
    dlg = QMessageBox(app)
    dlg.setWindowTitle("Importing file")
    file_path, _ = QFileDialog.getOpenFileName(app, "Choose a file", "", "FASTA Files (*.fasta *.fa);;All Files (*)")
    if file_path:
        try:
            SeqIO.parse(file_path, "fasta")
        except Exception as e:
            filename_label.setText(f"Choosen file: None")
            dlg.setText(f"Couldn't import a file:\n{e}")
            dlg.exec()
            return None
        else:
            filename_label.setText(f"Choosen file: {file_path}")
            return file_path
    else:
        filename_label.setText(f"Choosen file: None")
        dlg.setText("Importing file was canceled")
        dlg.exec()
        return None
        
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
    formatted_string = " ".join(text[i:i+125] for i in range(0, len(text), 125))
    return formatted_string

def save_file(app, text):
    dlg = QMessageBox(app)
    dlg.setWindowTitle("Saving file")
    file_path, _ = QFileDialog.getSaveFileName(app, "Save file", "", "Text Files (*.txt);;All Files (*)")
    if file_path:
        try:
            with open(file_path, "w", encoding="utf-8") as file:
                file.write(text)
        except Exception as e:
            dlg.setText(f"Couldn't save a file: {e}")
        else:
            dlg.setText("File saved successfully")
    else:
        dlg.setText("Saving file was canceled")
    dlg.exec()
