from PySide6.QtCore import Qt
from PySide6.QtWidgets import  QApplication, QMainWindow, QFileDialog, QLabel, QWidget, QScrollArea, QComboBox, QPushButton, QCheckBox, QSpacerItem, QSizePolicy, QMessageBox
from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout
from menus_dicts import dict_to_object_dict, starting_operations
from file_operations import *
from sequence_operations import *
from other_operations import *

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.file_operations = {
            "Get number of records in file": lambda: records_num(self.file_path),
            "Get length of each record": lambda: records_length(self.file_path),
            "Get the longest and the sortest record": lambda: find_longest_shortest(self.file_path)
        }

        self.dna_operations = {
            "Get sequence length": lambda: seq_length(self.file_path, self.sequence_id),
            "Get number of each nucleotide": lambda: nucleotides_count(self.file_path, self.sequence_id, self.sequence_type),
            "Get GC %": lambda: gc_content(self.file_path, self.sequence_id),
            "Get complementary sequence": lambda: get_complementary(self.file_path, self.sequence_id),
            "Get reverse complementary sequence": lambda: get_reverse_complementary(self.file_path, self.sequence_id),
            "Transcribe to mRNA": lambda: transcribe_to_mrna(self.file_path, self.sequence_id),
        }

        self.protein_operations = {
            "Get molecular weight in kDa": lambda: get_molecular_weight(self.file_path, self.sequence_id),
            "Get amino acids occurrence in number and %": lambda: get_amino_acids_occurrence(self.file_path, self.sequence_id),
            "Get protein isoelectric point": lambda: find_isoelectric_point(self.file_path, self.sequence_id),
            "Check protein stability": lambda: check_stability(self.file_path, self.sequence_id),
        }

        self.checkboxes_obj_list = []

        self.main_menu_objects = dict_to_object_dict(self, self.action_handler, "main")
        self.file_path = None
        self.sequence_type = ""
        self.sequence_id = ""
        self.additional_menu_state = ""

        self.setWindowTitle("FASTA File Perser")

        # Containers
        main_page_container = QWidget(self)
        main_page_container.setMinimumHeight(500)
        main_page_container.setMinimumWidth(1000)
        main_page_container.setStyleSheet("background-color: yellow;") ##

        self.additional_menu_container = QWidget(self)
        self.additional_menu_container.setMinimumWidth(400)
        self.additional_menu_container.setStyleSheet("background-color: lightblue;")

        menu_in_menu_cont = QWidget(self)
        menu_in_menu_cont.setMinimumWidth(300)
        menu_in_menu_cont.setStyleSheet("background-color: red;")

        # Layouts
        self.window_layout = QHBoxLayout()
        bottom_layout = QHBoxLayout()
        main_page_layout = QVBoxLayout(main_page_container)
        self.additional_menu_layout = QVBoxLayout(self.additional_menu_container)
        self.menu_in_menu_lo = QVBoxLayout(menu_in_menu_cont)

        # Label showing chosen file 
        self.filename_label = QLabel(f"Chosen file: {self.file_path}", self)
        self.filename_label.setStyleSheet("font-size: 15px")
        main_page_layout.addWidget(self.filename_label)

        # Scroll area for output
        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        # Output label nested in scroll area
        self.output_label = QLabel(self)
        self.output_label.setWordWrap(True)
        self.output_label.setStyleSheet("padding: 5px")
        self.output_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        scroll_area.setWidget(self.output_label)
        main_page_layout.addWidget(scroll_area)

        # Copy button
        copy_button = QPushButton("Copy output", self)
        copy_button.setFixedWidth(150)
        bottom_layout.addWidget(copy_button)

        # Overwrite output checkbox
        self.save_output_checkbox = QCheckBox("Don't overwrite outputs", self)
        self.save_output_checkbox.setChecked(False)
        bottom_layout.addWidget(self.save_output_checkbox)
        main_page_layout.addLayout(bottom_layout)

        self.window_layout.addWidget(main_page_container)

        # Additional menu
        self.title_label = QLabel("Avaliable operations", self)
        self.title_label.setStyleSheet("font-size: 20px; font-weight: bold;")

        self.additional_menu_layout.addWidget(self.title_label, alignment=Qt.AlignmentFlag.AlignCenter)
        spacer = QSpacerItem(0, 20, QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.additional_menu_layout.addSpacerItem(spacer)

        self.menu_in_menu_lo.setSpacing(10)

        self.create_menu("main", starting_operations)

        self.additional_menu_layout.addWidget(menu_in_menu_cont, alignment=Qt.AlignmentFlag.AlignCenter)

        self.window_layout.addWidget(self.additional_menu_container)
        self.window_layout.setAlignment(self.additional_menu_container, Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignRight)

        # Setting menu
        menu = self.menuBar()

        menu.addActions(list(self.main_menu_objects.values()))

        container = QWidget()
        container.setLayout(self.window_layout)
        self.setCentralWidget(container)


    def update_output(self, text):
        if self.save_output_checkbox.isChecked() and self.output_label.text() != "":
            separator = "----------------------------------------------------------------------------------------------------------"
            self.output_label.setText(f"{self.output_label.text()}\n{separator}\n{text}")
        else:
            self.output_label.setText(text)

    def action_handler(self):
        action = self.sender()
        action_id = action.data()

        if action_id == "import":
            self.file_path = check_file(self)
            if self.file_path:
                self.output_label.setText("")
                self.filename_label.setText(f"Choosen file: {self.file_path}")
                self.create_menu("main", starting_operations)
                self.sequence_type = ""
            else:
                self.output_label.setText("Importing file was canceled")
                self.filename_label.setText(f"Choosen file: None")
                self.create_menu("main", starting_operations)

        elif action_id == "file_menu":
            self.create_menu("file", self.file_operations)

        elif action_id == "sequence_menu":
            self.create_menu("sequence", self.dna_operations)


    def create_menu(self, menu_type, operations):
        if not self.file_path and menu_type != "main":
            dlg = QMessageBox(self)
            dlg.setWindowTitle("No file")
            dlg.setText("Please import file first")
            dlg.exec()
            return

        if menu_type == self.additional_menu_state:
            return

        # Delete content
        self.title_label.setText(f"{(menu_type.title())} Operations")
        self.delete_layout_content(0, self.menu_in_menu_lo)

        # Update with new one
        if menu_type == "main":
            for operation in operations:
                label = QLabel(operation, self)
                label.setWordWrap(True)
                label.setStyleSheet("font-size: 12px")
                self.menu_in_menu_lo.addWidget(label, alignment=Qt.AlignmentFlag.AlignLeft)

        elif menu_type == "sequence":
            self.sequence_combobox = QComboBox(self)
            self.sequence_combobox.addItem("Choose sequence from file")
            ids_list = records_ids(self.file_path)
            self.sequence_combobox.addItems(ids_list)
            self.sequence_combobox.currentIndexChanged.connect(self.sequence_changed)
            self.menu_in_menu_lo.addWidget(self.sequence_combobox)

        else:
            self.update_checkbox_list(operations)

        self.additional_menu_state = menu_type

    def sequence_changed(self, index):
            if self.sequence_combobox.itemText(0) == "Choose sequence from file" and index != 0:
                self.sequence_combobox.removeItem(0)
                index -= 1
            
            self.sequence_id = self.sequence_combobox.itemText(index)

            if self.sequence_type != check_seq_type(self.file_path, self.sequence_id):
                self.sequence_type = check_seq_type(self.file_path, self.sequence_id)
                self.delete_layout_content(1, self.menu_in_menu_lo)
                sequence_type_label = QLabel(f"Sequence type: {self.sequence_type}")
                self.menu_in_menu_lo.addWidget(sequence_type_label)
                if self.sequence_type == "DNA":
                    self.update_checkbox_list(self.dna_operations)
                elif self.sequence_type == "Protein":
                    self.update_checkbox_list(self.protein_operations)

    def delete_layout_content(self, index, layout):
        while layout.count():
            child = layout.takeAt(index)
            try: 
                if child.widget():
                    child.widget().deleteLater()
            except AttributeError:
                break

    def update_checkbox_list(self, operations):
        if len(self.checkboxes_obj_list) > 0:
                self.checkboxes_obj_list.clear()

        for name in operations.keys():
            checkbox = QCheckBox(name)
            self.menu_in_menu_lo.addWidget(checkbox)
            self.checkboxes_obj_list.append(checkbox)
        procees_buttoon = QPushButton("Proceed Operations")
        procees_buttoon.clicked.connect(lambda: self.process_selected_actions(operations))
        self.menu_in_menu_lo.addWidget(procees_buttoon)

    def process_selected_actions(self, operations):
        selected_operations = [
            cb.text() for cb in self.checkboxes_obj_list if cb.isChecked()
        ]

        if not selected_operations:
            self.output_label.setText("No operations selected.")
        else:
            result_data = [operations[op]() for op in selected_operations]
            self.output_label.setText(f"{'\n\n'.join(result_data)}")
            
app = QApplication([])
w = MainWindow()
w.show()
app.exec()
