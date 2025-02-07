from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QLabel, QHBoxLayout
from functionality import validate_answer

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Wordle Game")
        self.setFixedSize(510,650)
        self.setStyleSheet("background-color: #121213")

        self.word = "crown"

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

        self.answer_lineedit = QLineEdit()
        self.answer_lineedit.setPlaceholderText("Enter your guess")
        self.answer_lineedit.setStyleSheet("background-color: white; font-size: 15px")

        guess_button =  QPushButton("Guess")
        guess_button.setStyleSheet("background-color: white")
        guess_button.clicked.connect(self.guess_submitted)

        answer_layout.addWidget(self.answer_lineedit)
        answer_layout.addWidget(guess_button)

        self.warning_label = QLabel()
        self.warning_label.setStyleSheet("color: red; font-size: 15px")
        self.warning_label.setMaximumHeight(20)

        main_layout.addLayout(rows_layout)
        main_layout.addLayout(answer_layout)
        main_layout.addWidget(self.warning_label)
        main_layout.setAlignment(self.warning_label, Qt.AlignmentFlag.AlignHCenter)

        self.setLayout(main_layout)

    def guess_submitted(self):
        is_validated = validate_answer(self.answer_lineedit.text().strip(), self.warning_label)
        if is_validated:
            pass

app = QApplication([])
w = MainWindow()
w.show()
app.exec()
