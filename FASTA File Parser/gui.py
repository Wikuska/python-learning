from Bio import SeqIO
from PySide6.QtCore import QSize, Qt
from PySide6.QtWidgets import  QApplication, QMainWindow, QFileDialog, QLabel, QWidget, QScrollArea
from PySide6.QtWidgets import QVBoxLayout
from PySide6.QtGui import QAction

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

        self.output_label = QLabel(self)
        self.output_label.setWordWrap(True)
        self.output_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        self.scroll_area.setWidget(self.output_label)

        layout.addWidget(self.scroll_area)
        
        menu = self.menuBar()

        import_action = QAction("&Import File", self)
        import_action.triggered.connect(self.import_file)
        import_file = menu.addMenu("&Import File")
        import_file.addAction(import_action)

        records_number = QAction("&Show number of records", self)
        file_operations = menu.addMenu("&File Operations")
        file_operations.addAction(records_number)

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

app = QApplication([])
w = MainWindow()
w.show()
app.exec()
