from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit

class PlayerNameDialog(QDialog):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Enter players names")

        buttons = QDialogButtonBox.Ok
        self.buttonBox = QDialogButtonBox(buttons)
        self.buttonBox.accepted.connect(self.accept)
        
        p1_layout = QVBoxLayout()
        p2_layout = QVBoxLayout()
        horizontal_layout = QHBoxLayout()
        main_layout = QVBoxLayout()
        
        p1_label = QLabel("Player 1 name:")
        self.p1_input = QLineEdit()
        self.p1_input.setPlaceholderText("Player 1")
        p2_label = QLabel("Player 2 name:")
        self.p2_input = QLineEdit()
        self.p2_input.setPlaceholderText("Player 2")
        
        p1_layout.addWidget(p1_label)
        p1_layout.addWidget(self.p1_input)
        p2_layout.addWidget(p2_label)
        p2_layout.addWidget(self.p2_input)
        
        horizontal_layout.addLayout(p1_layout)
        horizontal_layout.addLayout(p2_layout)
        
        main_layout.addLayout(horizontal_layout)
        main_layout.addWidget(self.buttonBox)
        
        self.setLayout(main_layout)