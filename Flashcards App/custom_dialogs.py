from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QLabel, QLineEdit

class AddSetDialog(QDialog):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Create set")

        QBtn = (
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        )

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        layout = QVBoxLayout()
        self.main_label = QLabel("Name your set:", alignment = Qt.AlignmentFlag.AlignHCenter)
        self.set_name_input = QLineEdit(alignment=Qt.AlignmentFlag.AlignHCenter)

        layout.addWidget(self.main_label, alignment=Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(self.set_name_input)
        layout.addWidget(self.buttonBox, alignment = Qt.AlignmentFlag.AlignHCenter)
        self.setLayout(layout)
