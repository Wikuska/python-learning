from Bio import SeqIO
from PySide6.QtCore import QSize, Qt
from PySide6.QtWidgets import  QApplication, QMainWindow, QFileDialog, QLabel, QWidget, QScrollArea, QComboBox, QPushButton
from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout
from PySide6.QtGui import QAction
from operations import records_ids, check_seq_type, records_num, records_length, find_longest_shortest
from sequence_operations import seq_length, nucleotides_count, gc_content, get_complementary, get_reverse_complementary, transcribe_seq

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.file_name = ""
        self.sequence_id = ""
        self.sequence_type = ""

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
        copy_button.clicked.connect(self.copy_text)
        main_layout.addWidget(copy_button)
        
        # Setting menu
        menu = self.menuBar()

        # Import file menu action
        import_action = QAction("&Import File", self)
        import_action.triggered.connect(self.import_file)
        import_file = menu.addMenu("&Import File")
        import_file.addAction(import_action)

        # File operations menu actions
        records_number_action = QAction("&Get number of records", self)
        records_number_action.triggered.connect(self.records_number_act)
        records_length_action = QAction("&Get length of each record", self)
        records_length_action.triggered.connect(self.records_length_act)
        shortest_longest_record_action = QAction("&Get the longest and the sortest record", self)
        shortest_longest_record_action.triggered.connect(self.shortest_longest_record_act)

        # Setting file operations menu
        self.file_operations = menu.addMenu("&File Operations")
        self.file_operations.setEnabled(False)
        self.file_operations.addActions([records_number_action, records_length_action, shortest_longest_record_action])

        # Setting sequence opertions menu
        self.sequence_operations = menu.addMenu("&Sequence Operations")
        self.sequence_operations.setEnabled(False)

        # DNA menu actions
        sequence_length_action = QAction("&Get sequence length", self)
        sequence_length_action.triggered.connect(self.sequence_length_act)
        nucleotides_count_action = QAction("&Get number of each nucleotide", self)
        nucleotides_count_action.triggered.connect(self.nucleotides_count_act)
        gc_content_action = QAction("&Get GC %", self)
        gc_content_action.triggered.connect(self.gc_content_act)
        complementary_action = QAction("&Get complementary sequence", self)
        complementary_action.triggered.connect(self.complementary_act)
        reverse_complementary_action = QAction("&Get reverse complementary sequence", self)
        reverse_complementary_action.triggered.connect(self.reverse_complementary_act)
        transcribe_action = QAction("&Transcribe to mRNA", self)
        transcribe_action.triggered.connect(self.transcribe_act)

        # Setting DNA operations menu
        dna_operations = self.sequence_operations.addMenu("&DNA")
        dna_operations.addActions([sequence_length_action, nucleotides_count_action, gc_content_action, complementary_action, reverse_complementary_action, transcribe_action])

        # Adding all widgets to main window
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

    # Imort file menu actions
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

                ids_list = records_ids(self.file_name)
                for i in range(self.choose_seq_combobox.count() -1, 0, -1):
                    self.choose_seq_combobox.removeItem(i)
                self.choose_seq_combobox.addItems(ids_list)
                self.choose_seq_combobox.setEnabled(True)

    # File operations actions
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

    # DNA operations actions
    def sequence_length_act(self):
        length = seq_length(self.file_name, self.sequence_id)
        self.output_label.setText(f"Length of this sequence: {length}")

    def nucleotides_count_act(self):
        count_dict = nucleotides_count(self.file_name, self.sequence_id, self.sequence_type)
        self.output_label.setText("")
        for key, value in count_dict.items():
            self.output_label.setText(f"{self.output_label.text()} {key}: {value}\n")

    def gc_content_act(self):
        content = gc_content(self.file_name, self.sequence_id)
        self.output_label.setText(f"GC % of this sequence: {content}%")

    def complementary_act(self):
        complementary_seq = get_complementary(self.file_name, self.sequence_id)
        self.output_label.setText(f"Complementary sequence:\n{complementary_seq}")

    def reverse_complementary_act(self):
        reverse_comp_seq = get_reverse_complementary(self.file_name, self.sequence_id)
        self.output_label.setText(f"Reverse complementary sequence:\n{reverse_comp_seq}")

    def transcribe_act(self):
        mrna = transcribe_seq(self.file_name, self.sequence_id)
        self.output_label.setText(f"Template RNA:\n{mrna}")

    # Combobox text change handling   
    def sequence_changed(self, text):
        if text == "None":
            self.sequence_id = ""
            self.sequence_type = ""
            self.sequence_operations.setEnabled(False)
            self.seq_type_label.setText("Sequence type: None")
        else:
            self.sequence_id = text
            self.sequence_type = check_seq_type(self.file_name, self.sequence_id)
            self.seq_type_label.setText(f"Sequence type: {self.sequence_type}")
            self.sequence_operations.setEnabled(True)

    # Copy button action
    def copy_text(self):
        clipboard = QApplication.clipboard()
        clipboard.setText(self.output_label.text())


# Run app
app = QApplication([])
w = MainWindow()
w.show()
app.exec()
