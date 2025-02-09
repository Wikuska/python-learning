from PySide6.QtCore import Qt
from PySide6.QtWidgets import  QApplication, QLineEdit, QLabel, QWidget, QScrollArea, QComboBox, QPushButton, QCheckBox, QSpacerItem, QSizePolicy, QMessageBox, QStackedWidget
from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout
from PySide6.QtWidgets import QStyleFactory
from functionality import get_random_card, count_score

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Flashcard App")
        self.setFixedSize(880, 510)
        self.setStyleSheet("background-color: #00994C")

        self.dealer_ranks = []
        self.player_ranks = []

        self.dealer_score = 0
        self.player_score = 0

        main_layout = QVBoxLayout(self)
        score_layout = QHBoxLayout()
        score_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.dealer_layout = QHBoxLayout()
        self.dealer_layout.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft)
        self.player_layout = QHBoxLayout()
        self.player_layout.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft)
        buttons_layout = QHBoxLayout()

        self.dealer_score_label = QLabel(f"Dealer Score: 0")
        self.dealer_score_label.setStyleSheet("font-size: 15px;")
        self.player_score_label = QLabel(f"Your Score: 0")
        self.player_score_label.setStyleSheet("font-size: 15px;")

        self.set_cards()

        hit_button = QPushButton("HIT")
        hit_button.clicked.connect(lambda: self.draw_card(1))
        stand_button = QPushButton("STAND")

        score_layout.addWidget(self.dealer_score_label)
        score_layout.addWidget(self.player_score_label)

        buttons_layout.addWidget(hit_button)
        buttons_layout.addWidget(stand_button)

        main_layout.addLayout(score_layout, stretch=0)
        main_layout.addLayout(self.dealer_layout, stretch=1)
        main_layout.addSpacing(20)
        main_layout.addLayout(self.player_layout, stretch=1)
        main_layout.addLayout(buttons_layout, stretch=0)
        main_layout.setAlignment(buttons_layout, Qt.AlignmentFlag.AlignVCenter)

    def set_cards(self):
        get_random_card(self.dealer_layout, True)
        self.draw_card(0)
        self.draw_card(1)

    def draw_card(self, player):
        if player == 1:
            self.player_ranks.append(get_random_card(self.player_layout, False))
            self.player_score = count_score(self.player_ranks)
        else:
            self.dealer_ranks.append(get_random_card(self.dealer_layout, False))
            self.dealer_score = count_score(self.dealer_ranks)
        self.update_score()

    def update_score(self):
        self.dealer_score_label.setText(f"Dealer Score: {self.dealer_score}")
        self.player_score_label.setText(f"Your Score: {self.player_score}")


if __name__ == "__main__":
    app = QApplication([])
    app.setStyle(QStyleFactory.create("Fusion"))
    window = MainWindow()
    window.show()
    app.exec()