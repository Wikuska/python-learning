from PySide6.QtCore import QSize, Qt
from PySide6.QtWidgets import  QApplication, QMainWindow, QFileDialog, QLabel, QWidget, QScrollArea, QComboBox, QPushButton
from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout
from menus_dicts import dict_to_object_dict
from file_actions import *
from other_operations import *
from sequence_actions import *

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.file_path = ""
        self.sequence_id = ""
        self.sequence_type = ""

        self.import_menu_objects = dict_to_object_dict(self, self.action_handler, "import")
        self.file_menu_objects = dict_to_object_dict(self, self.action_handler, "file")
        self.dna_menu_objects = dict_to_object_dict(self, self.action_handler, "dna")

        # Setting window
        self.setWindowTitle("FASTA File Perser")
        self.setFixedSize(QSize(600, 300))

        # Layouts
        main_layout = QVBoxLayout()
        choose_seq_layout = QHBoxLayout()

        # Label showing chosen file 
        self.filename_label = QLabel("Chosen file: None", self)
        main_layout.addWidget(self.filename_label)

        # Label and combobox for choosing sequence
        self.choose_seq_label = QLabel("Chosen sequence: ", self)
        choose_seq_layout.addWidget(self.choose_seq_label)

        self.choose_seq_combobox = QComboBox(self)
        self.choose_seq_combobox.setFixedWidth(150)
        self.choose_seq_combobox.addItem("None")
        self.choose_seq_combobox.setEnabled(False)
        self.choose_seq_combobox.currentTextChanged.connect(self.sequence_changed)

        choose_seq_layout.addWidget(self.choose_seq_combobox)

        # Label for chosen sequence type
        self.seq_type_label = QLabel("Sequence type: None", self)
        choose_seq_layout.addWidget(self.seq_type_label)

        choose_seq_layout.addStretch()
        main_layout.addLayout(choose_seq_layout)

        # Scroll area for output
        self.scroll_area = QScrollArea(self)
        self.scroll_area.setWidgetResizable(True)
        # Output label nested in scroll area
        self.output_label = QLabel("Please import your file", self)
        self.output_label.setWordWrap(True)
        self.output_label.setStyleSheet("padding: 5px")
        self.output_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        self.scroll_area.setWidget(self.output_label)
        main_layout.addWidget(self.scroll_area)

        # Copy button
        copy_button = QPushButton("Copy output", self)
        copy_button.clicked.connect(self.copy_button)
        main_layout.addWidget(copy_button)
        
        # Setting menu
        menu = self.menuBar()

        self.import_operations = menu.addMenu("&Import")
        self.import_operations.addActions(list(self.import_menu_objects.values()))

        self.file_operations = menu.addMenu("&File Operations")
        self.file_operations.addActions(list(self.file_menu_objects.values()))
        self.file_operations.setEnabled(False)


        self.sequence_operations = menu.addMenu("&Sequence Operations")
        self.sequence_operations.setEnabled(False)

        self.dna_operations = self.sequence_operations.addMenu("&DNA")
        self.dna_operations.addActions(list(self.dna_menu_objects.values()))
        
        # Adding all widgets to main window
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)


    def action_handler(self):
        action = self.sender()
        action_id = action.data()

        if action_id == "import_file":
            file_path, _ = QFileDialog.getOpenFileName(self, "Choose a file", "", "FASTA Files (*.fasta *.fa);;All Files (*)")
            if file_path:
                try:
                    SeqIO.parse(file_path, "fasta")
                except Exception as e:
                    self.output_label.setText(f"An unexpected error occured: {e}")
                else:
                    self.file_path = file_path
                    self.filename_label.setText(f"Choosen file: {self.file_path}")
                    self.file_operations.setEnabled(True)

                    ids_list = records_ids(self.file_path)
                    for i in range(self.choose_seq_combobox.count() -1, 0, -1):
                        self.choose_seq_combobox.removeItem(i)
                    self.choose_seq_combobox.addItems(ids_list)
                    self.choose_seq_combobox.setEnabled(True)

        elif action_id == "records_number":
            self.output_label.setText(records_num(self.file_path))

        elif action_id == "records_length":
            self.output_label.setText(records_length(self.file_path))

        elif action_id == "shortest_longest_record":
            self.output_label.setText(find_longest_shortest(self.file_path))

        elif action_id == "sequence_length":
            self.output_label.setText(seq_length(self.file_path, self.sequence_id))

        elif action_id == "nucleotides_count":
            self.output_label.setText(nucleotides_count(self.file_path, self.sequence_id, self.sequence_type))

        elif action_id == "gc_content":
            self.output_label.setText(gc_content(self.file_path, self.sequence_id))

        elif action_id == "complementary":
            self.output_label.setText(get_complementary(self.file_path, self.sequence_id))

        elif action_id == "reverse_complementary":
            self.output_label.setText(get_reverse_complementary(self.file_path, self.sequence_id))

        elif action_id == "transcribe":
            self.output_label.setText(transcribe_to_mrna(self.file_path, self.sequence_id))

    def sequence_changed(self, text):
        if text == "None":
            self.sequence_id = ""
            self.sequence_type = ""
            self.sequence_operations.setEnabled(False)
            self.seq_type_label.setText("Sequence type: None")
        else:
            self.sequence_id = text
            self.sequence_type = check_seq_type(self.file_path, self.sequence_id)
            self.seq_type_label.setText(f"Sequence type: {self.sequence_type}")
            self.sequence_operations.setEnabled(True)

    def copy_button(self):
        clipboard = QApplication.clipboard()
        clipboard.setText(self.output_label.text())


# Run app
app = QApplication([])
w = MainWindow()
w.show()
app.exec()
