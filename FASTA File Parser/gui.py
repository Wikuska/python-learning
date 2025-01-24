from Bio import SeqIO
from PySide6.QtCore import QSize, Qt
from PySide6.QtWidgets import  QApplication, QMainWindow, QFileDialog, QLabel, QWidget, QScrollArea
from PySide6.QtWidgets import QVBoxLayout
from PySide6.QtGui import QAction
from operations import records_num, records_length, find_longest_shortest

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.file_name = ""

        self.setWindowTitle("FASTA File Perser")
        self.setFixedSize(QSize(600, 300))

        layout = QVBoxLayout()

        self.filename_label = QLabel("Chosen file: None", self)
        layout.addWidget(self.filename_label)

        self.scroll_area = QScrollArea(self)
        self.scroll_area.setWidgetResizable(True)

        self.output_label = QLabel("Please import your file", self)
        self.output_label.setWordWrap(True)
        self.output_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        self.scroll_area.setWidget(self.output_label)

        layout.addWidget(self.scroll_area)
        
        menu = self.menuBar()

        import_action = QAction("&Import File", self)
        import_action.triggered.connect(self.import_file)
        import_file = menu.addMenu("&Import File")
        import_file.addAction(import_action)

        records_number_action = QAction("&Get number of records", self)
        records_number_action.triggered.connect(self.records_number_act)
        records_length_action = QAction("&Get length of each record", self)
        records_length_action.triggered.connect(self.records_length_act)
        shortest_longest_record_action = QAction("&Get the longest and the sortest record", self)
        shortest_longest_record_action.triggered.connect(self.shortest_longest_record_act)

        self.file_operations = menu.addMenu("&File Operations")
        self.file_operations.setEnabled(False)
        self.file_operations.addActions([records_number_action, records_length_action, shortest_longest_record_action])

        sequence_operations = menu.addMenu("&Sequence Operations")

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def import_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Choose a file", "", "FASTA Files (*.fasta *.fa);;All Files (*)")
        if file_path:
            try:
                SeqIO.parse(file_path, "fasta")
            except Exception as e:
                self.output_label.setText(f"An unexpected error occured: {e}")
            else:
                self.file_name = file_path
                self.filename_label.setText(f"Choosen file: {self.file_name}")
                self.output_label.setText("File successfully uploaded")
                self.file_operations.setEnabled(True)

    def records_number_act(self):
        records = records_num(self.file_name)
        if records > 1:
            self.output_label.setText(f"There are {records} records in this file")
        elif records == 1:
            self.output_label.setText(f"There is only {records} record in this file")

    def records_length_act(self):
        length_dict = records_length(self.file_name)
        self.output_label.setText("Processing items: \n")
        for key, value in length_dict.items():
            self.output_label.setText(self.output_label.text() + f"{key}. {value['description']} - {value['length']}\n")
            QApplication.processEvents()

        self.output_label.setText(self.output_label.text() + "All items processed!")

    def shortest_longest_record_act(self):
        longest_shortest_dict = find_longest_shortest(self.file_name)
        if not longest_shortest_dict:
            self.output_label.setText("There is only 1 record in this file")
        else:
            self.output_label.setText(f"Longest:\nID: {longest_shortest_dict["longest"]["id"]}\nDescription: {longest_shortest_dict["longest"]["description"]}\nLength: {longest_shortest_dict["longest"]["length"]}")
            self.output_label.setText(self.output_label.text() + f"\n\nShortest:\nID: {longest_shortest_dict["shortest"]["id"]}\nDescription: {longest_shortest_dict["shortest"]["description"]}\nLength: {longest_shortest_dict["shortest"]["length"]}")
        

app = QApplication([])
w = MainWindow()
w.show()
app.exec()
