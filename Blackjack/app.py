from PySide6.QtCore import Qt
from PySide6.QtWidgets import  QApplication, QLineEdit, QLabel, QWidget, QScrollArea, QComboBox, QPushButton, QCheckBox, QSpacerItem, QSizePolicy, QMessageBox, QStackedWidget
from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout
from PySide6.QtWidgets import QStyleFactory

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Flashcard App")
        self.setFixedSize(880, 510)
        self.setStyleSheet("background-color: #00994C")

        self.dealer_score = 0
        self.player_score = 0

        main_layout = QVBoxLayout(self)
        score_layout = QHBoxLayout()
        dealer_layout = QHBoxLayout()
        player_layout = QHBoxLayout()
        buttons_layout = QHBoxLayout()

        dealer_score_label = QLabel(f"Dealer Score: {self.dealer_score}")
        dealer_score_label.setStyleSheet("font-size: 15px;")
        player_score_label = QLabel(f"Your Score: {self.player_score}")
        player_score_label.setStyleSheet("font-size: 15px;")

        example_dealer_card = QLabel("6 ♡", alignment = Qt.AlignmentFlag.AlignCenter)
        example_dealer_card.setFixedSize(100,170)
        example_dealer_card.setStyleSheet("font-size: 20px; color: red; background-color: white")
        example_dealer_card.setFrameShape(QLabel.Box)
        example_dealer_card.setLineWidth(3)

        example_player_card = QLabel("K ♧", alignment = Qt.AlignmentFlag.AlignCenter)
        example_player_card.setFixedSize(100,170)
        example_player_card.setStyleSheet("font-size: 20px; color: black; background-color: white")
        example_player_card.setFrameShape(QLabel.Box)
        example_player_card.setLineWidth(3)

        hit_button = QPushButton("HIT")
        stand_button = QPushButton("STAND")

        score_layout.addWidget(dealer_score_label)
        score_layout.addWidget(player_score_label)
        score_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)

        dealer_layout.addWidget(example_dealer_card, alignment = Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft)

        player_layout.addWidget(example_player_card, alignment = Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft)

        buttons_layout.addWidget(hit_button)
        buttons_layout.addWidget(stand_button)

        main_layout.addLayout(score_layout, stretch=0)
        main_layout.addLayout(dealer_layout, stretch=1)
        main_layout.addSpacing(20)
        main_layout.addLayout(player_layout, stretch=1)
        main_layout.addLayout(buttons_layout, stretch=0)
        main_layout.setAlignment(buttons_layout, Qt.AlignmentFlag.AlignVCenter)


if __name__ == "__main__":
    app = QApplication([])
    app.setStyle(QStyleFactory.create("Fusion"))
    window = MainWindow()
    window.show()
    app.exec()