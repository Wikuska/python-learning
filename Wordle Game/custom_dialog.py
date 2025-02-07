from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QLabel

class CustomDialog(QDialog):
    def __init__(self, is_won, word):
        super().__init__()
        
        if is_won:
            title = "Congratulations!"
            msg = "You guessed correct word!\nWanna play again?"
        else:
            title = "Try Again!"
            msg = f"Correct word was: {word}\nWanna play again?"

        self.setWindowTitle(title)

        QBtn = (
            QDialogButtonBox.Yes | QDialogButtonBox.No
        )

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        layout = QVBoxLayout()
        message = QLabel(msg, alignment = Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(message, alignment=Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(self.buttonBox, alignment = Qt.AlignmentFlag.AlignHCenter)
        self.setLayout(layout)
