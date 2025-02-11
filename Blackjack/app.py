from PySide6.QtCore import Qt
from PySide6.QtWidgets import  QApplication, QLineEdit, QLabel, QWidget, QScrollArea, QComboBox, QPushButton, QCheckBox, QSpacerItem, QSizePolicy, QMessageBox, QStackedWidget
from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout
from PySide6.QtWidgets import QStyleFactory
from functionality import get_random_card, count_score, get_hidden_card, clear_layout
from custom_dialog import CustomDialog

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.suits = ['♠', '♣', '♦', '♥']
        self.values = ['2', '3', '4', '5', '6', '7', '8', '9', '10', 'J', 'Q', 'K', 'A']

        self.setWindowTitle("Flashcard App")
        self.setFixedSize(880, 510)
        self.setStyleSheet("background-color: #00994C")

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

        self.hit_button = QPushButton("HIT")
        self.hit_button.clicked.connect(lambda: self.draw_card(1))
        stand_button = QPushButton("STAND")
        stand_button.clicked.connect(self.end_players_turn)

        score_layout.addWidget(self.dealer_score_label)
        score_layout.addWidget(self.player_score_label)

        buttons_layout.addWidget(self.hit_button)
        buttons_layout.addWidget(stand_button)

        self.set_up_game()

        main_layout.addLayout(score_layout, stretch=0)
        main_layout.addLayout(self.dealer_layout, stretch=1)
        main_layout.addSpacing(20)
        main_layout.addLayout(self.player_layout, stretch=1)
        main_layout.addLayout(buttons_layout, stretch=0)
        main_layout.setAlignment(buttons_layout, Qt.AlignmentFlag.AlignVCenter)

    def set_up_game(self):
        self.hit_button.setEnabled(True)

        self.deck = [f"{value} {suit}" for suit in self.suits for value in self.values]

        self.dealer_ranks = []
        self.player_ranks = []
        self.dealer_score = 0
        self.player_score = 0

        clear_layout(self.dealer_layout)
        clear_layout(self.player_layout)

        self.hidden_card = get_random_card(self.dealer_layout, True, self.deck)
        self.draw_card(0)
        self.draw_card(1)

    def draw_card(self, player):
        print(len(self.deck))
        if player == 1:
            self.player_ranks.append(get_random_card(self.player_layout, False, self.deck))
            self.player_score = count_score(self.player_ranks)
        else:
            self.dealer_ranks.append(get_random_card(self.dealer_layout, False, self.deck))
            self.dealer_score = count_score(self.dealer_ranks)
        self.update_score()

    def update_score(self):
        self.dealer_score_label.setText(f"Dealer Score: {self.dealer_score}")
        self.player_score_label.setText(f"Your Score: {self.player_score}")
        if self.player_score > 21:
            self.game_ends(0, "Your score is over 21!")
        elif self.dealer_score > 21:
            self.game_ends(1, "Dealer score is over 21!")


    def end_players_turn(self):
        self.dealer_ranks.append(get_hidden_card(self.hidden_card, self.deck))
        self.dealer_score = count_score(self.dealer_ranks)
        self.update_score()
        self.hit_button.setEnabled(False)

        while self.dealer_score < 17:
            self.draw_card(0)

        if self.dealer_score > self.player_score:
            self.game_ends(0, "Dealer score is higher!")
        elif self.dealer_score < self.player_score:
            self.game_ends(1, "Your score is higher!")
        elif self.dealer_score == self.player_score:
            self.game_ends(2, "Scores are the same!")

    def game_ends(self, is_won, msg):
        dlg = CustomDialog(is_won, msg)
        if dlg.exec():
            self.set_up_game()
        else:
            self.close()


if __name__ == "__main__":
    app = QApplication([])
    app.setStyle(QStyleFactory.create("Fusion"))
    window = MainWindow()
    window.show()
    app.exec()