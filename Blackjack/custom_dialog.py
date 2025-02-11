from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QLabel

class CustomDialog(QDialog):
    def __init__(self, is_won, msg):
        super().__init__()
        
        if is_won == 1:
            title = "You won!"
        elif is_won == 0:
            title = "You lost!"
        else:
            title = "It's a tie!"

        self.setWindowTitle(title)

        QBtn = (
            QDialogButtonBox.Yes | QDialogButtonBox.No
        )

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        layout = QVBoxLayout()
        message = QLabel(f"{msg}\nWanna play again?", alignment = Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(message, alignment=Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(self.buttonBox, alignment = Qt.AlignmentFlag.AlignHCenter)
        self.setLayout(layout)
