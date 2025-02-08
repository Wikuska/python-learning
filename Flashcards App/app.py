from PySide6.QtCore import Qt
from PySide6.QtWidgets import  QApplication, QMainWindow, QFrame, QLabel, QWidget, QScrollArea, QComboBox, QPushButton, QCheckBox, QSpacerItem, QSizePolicy, QMessageBox, QStackedWidget
from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout
import sys
from functionality import get_topics, get_topic_question, get_question_answer, set_question_known

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Flashcard App")
        self.setFixedSize(430, 550)

        self.current_question = None

        # Create container for multiple screens
        self.stack = QStackedWidget(self)

        # Create pages
        self.menu_screen = self.create_main_screen()
        self.flashcard_screen = self.create_flashcard_screen()

        # Add pages to the stack
        self.stack.addWidget(self.menu_screen)
        self.stack.addWidget(self.flashcard_screen)

        layout = QVBoxLayout(self)
        layout.addWidget(self.stack)


    def create_main_screen(self):
        screen = QWidget()
        layout = QVBoxLayout(screen)

        start_button = QPushButton("Start learning")
        modify_button = QPushButton("Modify Data")

        start_button.clicked.connect(lambda: self.stack.setCurrentWidget(self.flashcard_screen))

        layout.addWidget(start_button)
        layout.addWidget(modify_button)

        return screen

    def create_flashcard_screen(self):
        screen = QWidget()
        layout = QVBoxLayout(screen)
        buttons_layout = QHBoxLayout()

        self.flashcard_topics_combobox = QComboBox()
        self.flashcard_topics_combobox.addItem("Choose topic")
        topics_list = get_topics()
        if topics_list:
            self.flashcard_topics_combobox.addItems(topics_list)
        self.flashcard_topics_combobox.currentIndexChanged.connect(self.topic_changed)

        self.include_known_checkbox = QCheckBox("Include questions marked as known")

        self.flashcard_label = QLabel("Question", alignment = Qt.AlignmentFlag.AlignCenter)
        self.flashcard_label.setWordWrap(True)
        self.flashcard_label.setStyleSheet("font-size: 25px; background-color: #E6E6E6")

        see_answer_button = QPushButton("Show me answer")
        see_answer_button.clicked.connect(self.show_answer)

        learned_button = QPushButton("I know that one")
        learned_button.clicked.connect(self.answer_known)

        need_to_practice_button = QPushButton("Need to practice")
        need_to_practice_button.clicked.connect(self.answer_not_known)

        buttons_layout.addWidget(learned_button)
        buttons_layout.addWidget(see_answer_button)
        buttons_layout.addWidget(need_to_practice_button)

        layout.addWidget(self.flashcard_topics_combobox)
        layout.addWidget(self.include_known_checkbox)
        layout.addWidget(self.flashcard_label)
        layout.addLayout(buttons_layout)

        return screen
    
    def topic_changed(self, index):
        if self.flashcard_topics_combobox.itemText(0) == "Choose topic" and index != 0:
            self.flashcard_topics_combobox.removeItem(0)
            index -= 1

        self.show_question()

    def show_question(self):
        self.current_question = get_topic_question(self.flashcard_topics_combobox.currentText(), self.current_question, self.include_known_checkbox.isChecked())

        if not self.current_question:
            self.flashcard_label.setText("You know everything from this set")
            return

        self.flashcard_label.setText(self.current_question)

    def show_answer(self):
        answer = get_question_answer(self.current_question)
        if answer:
            self.flashcard_label.setText(answer)

    def answer_known(self):
        set_question_known(1, self.current_question)
        self.show_question()

    def answer_not_known(self):
        set_question_known(0, self.current_question)
        self.show_question()
    

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
