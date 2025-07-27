from PySide6.QtWidgets import QGridLayout, QPushButton, QVBoxLayout, QLabel
from PySide6.QtCore import QSize

def create_board(name, size=10):
    grid_layout = QGridLayout()
    player_layout = QVBoxLayout()
    buttons = []
    
    name_label = QLabel(name)
    player_layout.addWidget(name_label)
    
    for row in range(size):
        row_buttons = []
        for col in range(size):
            button = QPushButton("")
            button.setFixedSize(QSize(40, 40))
            button.setStyleSheet("background-color: lightblue; border: 1px solid black;")
            grid_layout.addWidget(button, row, col)
            row_buttons.append(button)
        buttons.append(row_buttons)
    player_layout.addLayout(grid_layout)
    
    return player_layout, buttons