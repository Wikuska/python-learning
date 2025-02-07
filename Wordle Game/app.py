from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QLabel, QHBoxLayout

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Wordle Game")
        self.setFixedSize(510,650)
        self.setStyleSheet("background-color: #121213")

        label_objets = []

        main_layout = QVBoxLayout()
        rows_layout = QVBoxLayout()
        answer_layout = QHBoxLayout()

        for i in range(6):
            row_layout = QHBoxLayout()
            for i in range(5):
                label = QLabel(alignment = Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
                label.setStyleSheet("font-size: 24px; font-weight: bold; background-color: #3a3a3c")
                label_objets.append(label)
                row_layout.addWidget(label)
            rows_layout.addLayout(row_layout)

        answer_lineedit = QLineEdit()
        answer_lineedit.setPlaceholderText("Enter your guess")
        answer_lineedit.setStyleSheet("background-color: white; font-size: 15px")

        guess_button =  QPushButton("Guess")
        guess_button.setStyleSheet("background-color: white")

        answer_layout.addWidget(answer_lineedit)
        answer_layout.addWidget(guess_button)

        main_layout.addLayout(rows_layout)
        main_layout.addLayout(answer_layout)

        self.setLayout(main_layout)


app = QApplication([])
w = MainWindow()
w.show()
app.exec()
