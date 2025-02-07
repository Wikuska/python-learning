from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QLabel, QHBoxLayout
from functionality import validate_answer, check_correct_letters, create_labels
from custom_dialog import CustomDialog

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Wordle Game")
        self.setFixedSize(510,650)
        self.setStyleSheet("background-color: #121213")

        self.word = "crown"

        self.main_layout = QVBoxLayout()
        answer_layout = QHBoxLayout()

        self.label_objects = create_labels(self.main_layout)

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

        self.main_layout.addLayout(answer_layout)
        self.main_layout.addWidget(self.warning_label)
        self.main_layout.setAlignment(self.warning_label, Qt.AlignmentFlag.AlignHCenter)

        self.setLayout(self.main_layout)

    def guess_submitted(self):
        guess = self.answer_lineedit.text().strip().lower()
        is_validated = validate_answer(guess, self.warning_label)
        if is_validated:
            self.answer_lineedit.clear()
            game_won = check_correct_letters(guess, self.word, self.label_objects)

            if game_won or not game_won and len(self.label_objects) < 5:
                dlg = CustomDialog(game_won, self.word)
                if dlg.exec():
                    self.label_objects = create_labels(self.main_layout)
                else:
                    self.close()             

app = QApplication([])
w = MainWindow()
w.show()
app.exec()
